#!python
#cython: embedsignature=True
cimport numpy as np
from libc.stdio cimport FILE
from libc.stdlib cimport free


cdef extern from "common.h":
    ctypedef np.int32_t boolean
    ctypedef np.uint16_t bits16
    ctypedef np.uint32_t bits32
    ctypedef np.uint64_t bits64
    cdef struct slName:
        slName *next
        char name[1]

    int sameString(char *a, char *b)
    FILE *mustOpen(char *fileName, char *mode)
    void carefulClose(FILE **pFile)
    void slFreeList(void *listPt)
    void freeMem(void *pt)
    void AllocArray(void *pt, size_t size)


cdef extern from "localmem.h":
    cdef struct lm:
        size_t blockSize
        size_t allignMask
        size_t allignAdd

    lm *lmInit(int blockSize)
    void lmCleanup(lm **pLm)


cdef extern from "bPlusTree.h":
    cdef struct bptFile:
        bptFile *next
        udcFile *udc
        char *fileName
        bits32 blockSize
        bits32 keySize
        bits32 valSize
        bits64 itemCount
        boolean isSwapped
        bits64 rootOffset


cdef extern from "crTree.h":
    cdef struct cirTreeFile:
        cirTreeFile *next
        udcFile *udc
        char *fileName
        boolean isSwapped
        bits64 rootOffset
        bits32 blockSize
        bits64 itemCount
        bits32 startChromIx
        bits32 startBase
        bits32 endChromIx
        bits32 endBase
        bits64 fileSize
        bits32 itemsPerSlot


cdef extern from "udc.h":
    cdef struct udcFile:
        pass

    char *udcDefaultDir()
    udcFile *udcFileMayOpen(char *url, char *cacheDir)
    udcFile *udcFileOpen(char *url, char *cacheDir)
    void udcFileClose(udcFile **pFile)


cdef extern from "sig.h":
    bits32 bigWigSig
    bits32 bigBedSig


cdef extern from "bbiFile.h":
    cdef struct bbiFile:
        bbiFile *next
        char *fileName
        udcFile *udc
        bptFile *chromBpt
        cirTreeFile *unzoomedCir
        bbiZoomLevel *levelList
        boolean isSwapped
        bits32 typeSig
        bits16 version
        bits16 zoomLevels
        bits64 chromTreeOffset
        bits64 unzoomedDataOffset
        bits64 unzoomedIndexOffset
        bits16 fieldCount
        bits16 definedFieldCount
        bits64 asOffset
        bits64 totalSummaryOffset
        bits32 uncompressBufSize
        bits64 extensionOffset
        bits16 extensionSize
        bits16 extraIndexCount
        bits64 extraIndexListOffset

    cdef struct bbiChromInfo:
        bbiChromInfo *next
        char *name
        bits32 id, size

    cdef struct bbiInterval:
        bbiInterval *next
        bits32 start, end
        double val

    cdef struct bbiZoomLevel:
        bbiZoomLevel *next
        bits32 reductionLevel
        bits32 reserved
        bits64 dataOffset
        bits64 indexOffset

    cdef enum bbiSummaryType:
        bbiSumMean = 0
        bbiSumMax = 1
        bbiSumMin = 2
        bbiSumCoverage = 3
        bbiSumStandardDeviation = 4
        bbiSumSum = 5

    cdef struct bbiSummary:
        bbiSummary *next
        bits32 chromId
        bits32 start,end
        bits32 validCount
        float minVal
        float maxVal
        float sumData
        float sumSquares
        bits64 fileOffset

    cdef struct bbiSummaryElement:
        bits64 validCount
        double minVal
        double maxVal
        double sumData
        double sumSquares

    ctypedef bbiInterval *(*BbiFetchIntervals)(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        lm *lm)

    bbiFile *bbiFileOpen(char *fileName, bits32 sig, char *typeName)
    void bbiFileClose(bbiFile **pBwf)
    boolean bbiFileCheckSigs(char *fileName, bits32 sig, char *typeName)
    bits32 bbiChromSize(bbiFile *bbi, char *chrom)
    bbiChromInfo *bbiChromList(bbiFile *bbi)
    void bbiChromInfoFreeList(bbiChromInfo **pList)
    bbiZoomLevel *bbiBestZoom(
        bbiZoomLevel *levelList,
        int desiredReduction)
    bits32 bbiIntervalSlice(
        bbiFile *bbi,
        bits32 baseStart,
        bits32 baseEnd,
        bbiInterval *intervalList,
        bbiSummaryElement *el)
    bits32 bbiSummarySlice(
        bbiFile *bbi,
        bits32 baseStart,
        bits32 baseEnd,
        bbiSummary *sumList,
        bbiSummaryElement *el)
    bbiSummary *bbiSummariesInRegion(
        bbiZoomLevel *zoom,
        bbiFile *bbi,
        int chromId,
        bits32 start,
        bits32 end)
    boolean bbiSummaryArrayFromZoom(
        bbiZoomLevel *zoom,
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        int summarySize,
        bbiSummaryElement *summary)
    boolean bbiSummaryArrayFromFull(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        BbiFetchIntervals fetchIntervals,
        int summarySize,
        bbiSummaryElement *summary)
    boolean bbiSummaryArray(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        BbiFetchIntervals fetchIntervals,
        bbiSummaryType summaryType,
        int summarySize,
        double *summaryValues)
    boolean bbiSummaryArrayExtended(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        BbiFetchIntervals fetchIntervals,
        int summarySize,
        bbiSummaryElement *summary)
    bbiSummaryElement bbiTotalSummary(
        bbiFile *bbi)

    # Functions that were declared static (private) in bbiRead.c
    # Source was modified to expose them here via bbiFile.h
    int bbiChromId(bbiFile *bbi, char *chrom)
    bits32 bbiIntervalSlice(
        bbiFile *bbi,
        bits32 baseStart,
        bits32 baseEnd,
        bbiInterval *intervalList,
        bbiSummaryElement *el)
    bits32 bbiSummarySlice(
        bbiFile *bbi,
        bits32 baseStart,
        bits32 baseEnd,
        bbiSummary *sumList,
        bbiSummaryElement *el)

    # Functions that were declared static (private) in bigBed.c
    # Source was modified to expose them here via bbiFile.h
    bbiInterval *bigBedCoverageIntervals(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        lm *lm)


cdef extern from "bigWig.h":
    boolean isBigWig(char *fileName)
    bbiFile *bigWigFileOpen(char *fileName)
    bbiInterval *bigWigIntervalQuery(
        bbiFile *bwf,
        char *chrom,
        bits32 start,
        bits32 end,
        lm *lm)
    double bigWigSingleSummary(
        bbiFile *bwf,
        char *chrom,
        int start,
        int end,
        bbiSummaryType summaryType,
        double defaultVal)
    boolean bigWigSummaryArray(
        bbiFile *bwf,
        char *chrom,
        bits32 start,
        bits32 end,
        bbiSummaryType summaryType,
        int summarySize,
        double *summaryValues)
    boolean bigWigSummaryArrayExtended(
        bbiFile *bwf,
        char *chrom,
        bits32 start,
        bits32 end,
        int summarySize,
        bbiSummaryElement *summary)


cdef extern from "bigBed.h":
    cdef struct bigBedInterval:
        bigBedInterval *next
        bits32 start, end
        char *rest
        bits32 chromId

    bbiFile *bigBedFileOpen(char *fileName)
    #void bigBedFileClose(bbiFile **pBwf)
    bits64 bigBedItemCount(bbiFile *bbi)
    bigBedInterval *bigBedIntervalQuery(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        int maxItems,
        lm *lm)
    int bigBedIntervalToRow(
        bigBedInterval *interval,
        char *chrom,
        char *startBuf,
        char *endBuf,
        char **row,
        int rowSize)
    boolean bigBedSummaryArray(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        bbiSummaryType summaryType,
        int summarySize,
        double *summaryValues)
    boolean bigBedSummaryArrayExtended(
        bbiFile *bbi,
        char *chrom,
        bits32 start,
        bits32 end,
        int summarySize,
        bbiSummaryElement *summary)
    char *bbiCachedChromLookup(
        bbiFile *bbi,
        int chromId,
        int lastChromId,
        char *chromBuf,
        int chromBufSize)
    char *bigBedAutoSqlText(bbiFile *bbi)
    asObject *bigBedAs(bbiFile *bbi)
    asObject *bigBedAsOrDefault(bbiFile *bbi)
    asObject *bigBedFileAsObjOrDefault(char *fileName)


cdef extern from "asParse.h":
    cdef enum asTypes:
        t_double, t_float, t_char, t_int, t_uint, t_short, \
        t_ushort, t_byte, t_ubyte, t_off, t_string, t_lstring, t_object, \
        t_simple, t_enum, t_set

    cdef struct asTypeInfo:
        char *name
        bint isUnsigned
        bint stringy
        char *sqlName
        char *cName
        char *listyName
        char *nummyName
        char *outFormat
        char *djangoName

    cdef struct asIndex:
        asIndex *next
        char *type
        int size

    cdef struct asColumn:
        asColumn *next
        char *name
        char *comment
        asTypeInfo *lowType
        char *obName
        asObject *obType
        int fixedSize
        char *linkedSizeName
        asColumn *linkedSize
        bint isSizeLink
        bint isList
        bint isArray
        bint autoIncrement
        slName *values
        asIndex *index

    cdef struct asObject:
        asObject *next
        char *name
        char *comment
        asColumn *columnList
        bint isTable
        bint isSimple
        int bla

    asObject *asParseText(char *text)
    void asObjectFree(asObject **pAs)


cdef extern from "dystring.h":

    cdef struct dyString:
        dyString *next
        char *string
        int bufSize
        int stringSize


cdef extern from 'setjmp.h':
    struct __jmp_buf_tag:
        pass
    ctypedef __jmp_buf_tag jmp_buf
    int setjmp (jmp_buf __env)
    int longjmp (jmp_buf __env, int val)


cdef extern from "errCatch.h":

    cdef struct errCatch:
        errCatch *next
        jmp_buf jmpBuf
        dyString *message
        boolean gotError
        boolean gotWarning

    errCatch *errCatchNew()
    boolean errCatchStart(errCatch *errCatch)
    boolean errCatchPushHandlers(errCatch *errCatch)
    void errCatchEnd(errCatch *errCatch)
    void errCatchFree(errCatch *errCatch)


cdef extern from "basicBed.h":

    char *bedAsDef(int bedFieldCount, int totalFieldCount)
