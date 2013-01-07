/***************************************************************************
 * obspy-readbuffer.c:
 *
 * Reads a memory buffer to a MSTraceList structure and parses Selections.
 *
 * Parts are copied from tracelist.c and unpack.c from libmseed by Chad
 * Trabant
 *
 * There were two reasons why we chose to write our own C code using lower
 * level libmseed functions:
 * 
 * 1. ms_readtracelist() and ms_readtraces() only work with files and not
 * memory buffers. ObsPy only passes memory buffers to libmseed - we do this
 * even for files so we that we can have a unified handling whether the
 * MiniSEED data comes from the filesystem, internet or some other source. The
 * two routines listed above could probably be adopted quite easily to also do
 * this. <update> This issue is now resolved by the use of
 * msr_parse_selection </update>.
 * 
 * 2. realloc() is really slow on some Windows platforms [1] to the point
 * of becoming unusable after repeated calls. I think both the aforementioned
 * functions use it to continuously grow the data array when combining
 * records. This works just fine on Unix so it is a non-issue in most cases
 * but we would like to be able to support Windows. In our C code we first
 * read all records, then sort them, and then allocate a chunk of memory for
 * each continuous trace. This is also fast on Windows at the expense of
 * temporarily doubling the memory consumption.
 *
 * [1] http://blog.kowalczyk.info/article/2be/realloc-on-Windows-vs-Linux.html
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "libmseed/libmseed.h"


// Linkable container of MSRecords
typedef struct LinkedRecordList_s {
    struct MSRecord_s      *record;       // This record
    struct LinkedRecordList_s  *next;     // The next record container
}
LinkedRecordList;

// Container for a continuous linked list of records.
typedef struct ContinuousSegment_s {
    hptime_t starttime;                     // Time of the first sample
    hptime_t endtime;                       // Time of the last sample
    double samprate;                        // Sample rate
    char sampletype;                        // Sampletype
    hptime_t hpdelta;                       // High precission sample period
    int64_t samplecnt;                         // Total sample count
    /* Timing quality is a vendor specific value from 0 to 100% of maximum
     * accuracy, taking into account both clock quality and data flags. */
    uint8_t timing_qual;
    /* type of calibration available, BLK 300 = 1, BLK 310 = 2, BLK 320 = 3
     * BLK 390 = 4, BLK 395 = -2 */
    int8_t calibration_type;
    uint32_t blkt_buffer_len;
    uint8_t *blkt_buffer;
    void *datasamples;                      // Actual data samples
    struct LinkedRecordList_s *firstRecord; // First item
    struct LinkedRecordList_s *lastRecord;  // Last item
    struct ContinuousSegment_s *next;       // Next segment
    struct ContinuousSegment_s *previous;   // Previous segment
}
ContinuousSegment;

// A container for continuous segments with the same id
typedef struct LinkedIDList_s {
    char network[11];                         // Network designation, NULL terminated
    char station[11];                         // Station designation, NULL terminated
    char location[11];                        // Location designation, NULL terminated
    char channel[11];                         // Channel designation, NULL terminated
    char dataquality;                         // Data quality indicator */
    struct ContinuousSegment_s *firstSegment; // Pointer to first of list of segments
    struct ContinuousSegment_s *lastSegment;  // Pointer to last of list of segments
    struct LinkedIDList_s *next;              // Pointer to next id
    struct LinkedIDList_s *previous;          // Pointer to previous id
}
LinkedIDList;

// Init function for the LinkedIDList
static LinkedIDList *
lil_init(void)
{
    // Allocate 0 initialized memory.
    LinkedIDList *lil = (LinkedIDList *) malloc (sizeof(LinkedIDList));
    if ( lil == NULL ) {
        ms_log (2, "lil_init(): Cannot allocate memory\n");
        return NULL;
    }
    memset (lil, 0, sizeof (LinkedIDList));
    return lil;
}

// Init function for the LinkedRecordList
static LinkedRecordList *
lrl_init (void)
{
    // Allocate 0 initialized memory.
    LinkedRecordList *lrl = (LinkedRecordList *) malloc (sizeof(LinkedRecordList));
    if ( lrl == NULL ) {
        ms_log (2, "lrl_init(): Cannot allocate memory\n");
        return NULL;
    }
    memset (lrl, 0, sizeof (LinkedRecordList));
    return lrl;
}

// Init a Segment with a linked record list.
static ContinuousSegment *
seg_init(void)
{
    ContinuousSegment *seg = (ContinuousSegment *) malloc (sizeof(ContinuousSegment));
    if ( seg == NULL ) {
        ms_log (2, "seg_init(): Cannot allocate memory\n");
        return NULL;
    }
    memset (seg, 0, sizeof (ContinuousSegment));
    return seg;
}

// Frees a LinkedRecordList. The given Record is assumed to be the head of the
// list.
static void
lrl_free(LinkedRecordList * lrl)
{
    LinkedRecordList * next;
    while ( lrl != NULL) {
        next = lrl->next;
        msr_free(&lrl->record);
        free(lrl);
        if (next == NULL) {
            break;
        }
        lrl = next;
    }
    lrl = NULL;
}

// Frees a ContinuousSegment and all structures associated with it.
// The given segment is supposed to be the head of the linked list.
static void
seg_free(ContinuousSegment * seg)
{
    ContinuousSegment * next;
    while (seg != NULL) {
        next = seg->next;
        // free(seg->datasamples);
        if (seg->firstRecord != NULL) {
            lrl_free(seg->firstRecord);
        }
        free(seg);
        if (next == NULL) {
            break;
        }
        seg = next;
    }
    seg = NULL;
}

// Free a LinkedIDList and all structures associated with it.
void
lil_free(LinkedIDList * lil)
{
    LinkedIDList * next;
    while ( lil != NULL) {
        next = lil->next;
        if (lil->firstSegment != NULL) {
            seg_free(lil->firstSegment);
        }
        free(lil);
        if (next == NULL) {
            break;
        }
        lil = next;
    }
    lil = NULL;
}

static void
copySegmentData(ContinuousSegment * const contseg,
    const flag unpack_data, long (* const allocData) (int, char)) {
    int size;
    long offset;
    LinkedRecordList * reclst = NULL;

    if (contseg != NULL ) {
        if (contseg->datasamples) {
            free(contseg->datasamples);
        }
        // Allocate data via a callback function.
        if (unpack_data != 0) {
            contseg->datasamples = (void *) allocData(contseg->samplecnt,
                    contseg->sampletype);
        }

        // Loop over all records, write the data to the buffer and free the msr structures.
        reclst = contseg->firstRecord;
        offset = (long) (contseg->datasamples);
        while (reclst != NULL ) {
            size = reclst->record->samplecnt
                    * ms_samplesize(reclst->record->sampletype);
            memcpy((void *) offset, reclst->record->datasamples,
                    size);
            // Free the record.
            msr_free(&(reclst->record));
            // Increase the data_offset and the record.
            offset += (long) size;
            reclst = reclst->next;
        }
    }
}

typedef struct bfield_s {
    int32_t blkt_name;
    int32_t offset;
    int32_t size;
} binfield;

// Function that reads from a MiniSEED binary file from a char buffer and
// returns a LinkedIDList.
LinkedIDList *
readMSEEDBuffer (char *mseed, const int buflen, Selections *selections,
        const flag unpack_data, const int reclen, const flag verbose,
        const flag details, long (* const allocData) (int, char),
        const binfield * bfield, const int bfieldlen)
{
    // current offset of mseed char pointer
    int64_t offset = 0;
    /* TODO: change definition of readMSEEDBuffer */
    int64_t buflen64 = buflen;

    // the timing_qual of BLK 1001
    uint8_t timing_qual = 0xFF;

    // the calibration type, availability of BLK 300, 310, 320, 390, 395
    int8_t calibration_type = -1;

    // Init all the pointers to NULL. Most compilers should do this anyway.
    LinkedIDList * idListHead = NULL;
    LinkedIDList * idListCurrent = NULL;
    LinkedIDList * idListLast = NULL;
    MSRecord *msr = NULL;
    ContinuousSegment * segmentCurrent = NULL;
    hptime_t lastgap = 0;
    hptime_t hptimetol = 0;
    hptime_t nhptimetol = 0;
    LinkedRecordList *recordCurrent = NULL;
    int record_count = 0;
    uint32_t blkt_buffer_len = 0;
    uint8_t *blkt_buffer;
    uint32_t i;


    for (i = 0; i < bfieldlen; ++i) {
        blkt_buffer_len += bfield[i].size;
    }
    blkt_buffer = (uint8_t *) calloc(blkt_buffer_len, sizeof(uint8_t));

    /* Loop over all selected records in recbuf */
    while (offset < buflen64) {
        msr = msr_init(NULL);
        if (msr_parse_selection(mseed, buflen64, &offset, &msr, reclen,
                selections, unpack_data, verbose)) {
            msr_free(&msr);
            /* Only print error if offset is still within buffer length */
            if (verbose && offset < buflen64)
                ms_log(2, "Error parsing record at offset %lld0",
                        (long long int) offset);
            continue;
        }
        /* Successfully found and parsed record */
        record_count++;
        /* Increment offset in buffer for subsequent call to msr_parse_selection() */
        offset += msr->reclen;
        /* Do something with the record */
        recordCurrent = lrl_init ();
        // Append to linked record list if one exists.
        recordCurrent->record = msr;
        recordCurrent->next = NULL;

        // Check if the ID of the record is already available and if not create a
        // new one.
        // Start with the last id as it is most likely to be the correct one.
        idListCurrent = idListLast;
        while (idListCurrent != NULL) {
            if (strcmp(idListCurrent->network, recordCurrent->record->network) == 0 &&
                    strcmp(idListCurrent->station, recordCurrent->record->station) == 0 &&
                    strcmp(idListCurrent->location, recordCurrent->record->location) == 0 &&
                    strcmp(idListCurrent->channel, recordCurrent->record->channel) == 0 &&
                    idListCurrent->dataquality == recordCurrent->record->dataquality) {
                break;
            }
            else {
                idListCurrent = idListCurrent->previous;
            }
        }

        // Create a new id list if one is needed.
        if (idListCurrent == NULL) {
            idListCurrent = lil_init();
            idListCurrent->previous = idListLast;
            if (idListLast != NULL) {
                idListLast->next = idListCurrent;
            }
            idListLast = idListCurrent;
            if (idListHead == NULL) {
                idListHead = idListCurrent;
            }

            // Set the IdList attributes.
            strcpy(idListCurrent->network, recordCurrent->record->network);
            strcpy(idListCurrent->station, recordCurrent->record->station);
            strcpy(idListCurrent->location, recordCurrent->record->location);
            strcpy(idListCurrent->channel, recordCurrent->record->channel);
            idListCurrent->dataquality = recordCurrent->record->dataquality;
        }


        // Now check if the current record fits exactly to the end of the last
        // segment of the current id. If not create a new segment. Therefore
        // if records with the same id are in wrong order a new segment will be
        // created. This is on purpose.
        segmentCurrent = idListCurrent->lastSegment;
        if (segmentCurrent != NULL) {
            hptimetol = (hptime_t) (0.5 * segmentCurrent->hpdelta);
            nhptimetol = ( hptimetol ) ? -hptimetol : 0;
            lastgap = recordCurrent->record->starttime - segmentCurrent->endtime - segmentCurrent->hpdelta;
        }
        if ((details == 1) || (blkt_buffer_len >= 1)) {
            /* extract information on calibration BLKs */
            calibration_type = -1;
            if (recordCurrent->record->blkts) {
                BlktLink *cur_blkt = recordCurrent->record->blkts;
                while (cur_blkt) {
                    int step = 0;
                    for (i = 0; i < bfieldlen; ++i) {
                        if (cur_blkt->blkt_type == bfield[i].blkt_name) {
#if 0
                            printf("name %d, offset %d, size %d\n", bfield[i].blkt_name, bfield[i].offset, bfield[i].size);
#endif
                            memcpy((void *)(blkt_buffer + step), (void *)(cur_blkt->blktdata + bfield[i].offset), bfield[i].size);
                        }
                        step += bfield[i].size;
                    }
                    /* old way */
#if 0
                    printf("%d\n", cur_blkt->blkt_type);
#endif
                    switch (cur_blkt->blkt_type) {
                    case 300:
                        calibration_type = 1;
                        break;
                    case 310:
                        calibration_type = 2;
                        break;
                    case 320:
                        calibration_type = 3;
                        break;
                    case 390:
                        calibration_type = 4;
                        break;
                    case 395:
                        calibration_type = -2;
                        break;
                    default:
                        break;
                    }
                    cur_blkt = cur_blkt->next;
                }
            }
            /* extract information based on timing quality */
            timing_qual = 0xFF;
            if (recordCurrent->record->Blkt1001 != 0) {
                timing_qual = recordCurrent->record->Blkt1001->timing_qual;
            }
        }
        if ( segmentCurrent != NULL &&
             segmentCurrent->sampletype == recordCurrent->record->sampletype &&
             // Test the default sample rate tolerance: abs(1-sr1/sr2) < 0.0001
             MS_ISRATETOLERABLE (segmentCurrent->samprate, recordCurrent->record->samprate) &&
             // Check if the times are within the time tolerance
             lastgap <= hptimetol && lastgap >= nhptimetol &&
             segmentCurrent->timing_qual == timing_qual &&
             segmentCurrent->calibration_type == calibration_type &&
             memcmp(segmentCurrent->blkt_buffer, blkt_buffer, blkt_buffer_len) == 0) {
            segmentCurrent->lastRecord = segmentCurrent->lastRecord->next = recordCurrent;
            segmentCurrent->samplecnt += recordCurrent->record->samplecnt;
            segmentCurrent->endtime = msr_endtime(recordCurrent->record);
        }
        // Otherwise create a new segment and add the current record.
        else {
            // the last contiguous segment of the current can now be copied and
            // the corresponding records can be freed already
            copySegmentData(idListCurrent->lastSegment, unpack_data, allocData);
            segmentCurrent = seg_init();
            segmentCurrent->blkt_buffer = calloc(blkt_buffer_len, sizeof(uint8_t));
            memcpy(segmentCurrent->blkt_buffer, blkt_buffer, blkt_buffer_len);
            segmentCurrent->blkt_buffer_len = blkt_buffer_len;
            segmentCurrent->previous = idListCurrent->lastSegment;
            if (idListCurrent->lastSegment != NULL) {
                idListCurrent->lastSegment->next = segmentCurrent;
            }
            else {
                idListCurrent->firstSegment = segmentCurrent;
            }
            idListCurrent->lastSegment = segmentCurrent;

            segmentCurrent->starttime = recordCurrent->record->starttime;
            segmentCurrent->endtime = msr_endtime(recordCurrent->record);
            segmentCurrent->samprate = recordCurrent->record->samprate;
            segmentCurrent->sampletype = recordCurrent->record->sampletype;
            segmentCurrent->samplecnt = recordCurrent->record->samplecnt;
            // Calculate high-precision sample period
            segmentCurrent->hpdelta = (hptime_t) (( recordCurrent->record->samprate ) ?
                           (HPTMODULUS / recordCurrent->record->samprate) : 0.0);
            segmentCurrent->timing_qual = timing_qual;
            segmentCurrent->calibration_type = calibration_type;
            segmentCurrent->firstRecord = segmentCurrent->lastRecord = recordCurrent;
        }
    }
    // Return empty id list if no records could be found.
    if (record_count == 0) {
        idListHead = lil_init();
        return idListHead;
    }
    // copy data of all remaining segments
    idListCurrent = idListHead;
    while (idListCurrent != NULL) {
        copySegmentData(idListCurrent->lastSegment, unpack_data, allocData);
        idListCurrent = idListCurrent->next;
    }

    free(blkt_buffer);
    return idListHead;
}
