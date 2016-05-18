#!python
#cython: embedsignature=True
from __future__ import division, print_function
import numpy as np
import cython

cimport numpy as np
from libc.stdio cimport FILE


cdef extern from "common.h":
    ctypedef np.int32_t boolean
    ctypedef np.uint16_t bits16
    ctypedef np.uint32_t bits32
    ctypedef np.uint64_t bits64

    int sameString(char *a, char *b)
    FILE *mustOpen(char *fileName, char *mode)
    void carefulClose(FILE **pFile)


cdef extern from "localmem.h":
    cdef struct lm:
        size_t blockSize
        size_t allignMask
        size_t allignAdd

    lm *lmInit(int blockSize)
    void lmCleanup(lm **pLm)


cdef extern from "udc.h":
    cdef struct udcFile:
        pass


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


cdef extern from "bbiFile.h":
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

    cdef struct bbiChromInfo:
        bbiChromInfo *next
        char *name
        bits32 id, size

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
    
    bbiChromInfo *bbiChromList(bbiFile *bbi)
    void bbiFileClose(bbiFile **pBwf)
    void bbiChromInfoFreeList(bbiChromInfo **pList)
    

cdef extern from "bigWig.h":
    bbiFile *bigWigFileOpen(char *fileName)
    bbiInterval *bigWigIntervalQuery(bbiFile *bwf, char *chrom, bits32 start, bits32 end, lm *lm)



cdef inline _fetch(np.ndarray[np.double_t, ndim=1] out, 
            bbiFile *bwf, char *chromName, int start, int end,
            int refStart, int chromSize, double missing, double oob):
    cdef bbiInterval *intervalList
    cdef bbiInterval *interval
    cdef lm *lm
    cdef boolean firstTime
    cdef int saveStart, prevEnd
    cdef double saveVal

    lm = lmInit(0)
    intervalList = interval = bigWigIntervalQuery(bwf, chromName, start, end, lm)
    firstTime = True
    saveStart = -1
    prevEnd = -1
    saveVal = -1.0
    while interval != NULL:
        if firstTime:
            saveStart = interval.start
            saveVal = interval.val
            firstTime = False
        elif not ( (interval.start == prevEnd) and (interval.val == saveVal) ):
            out[saveStart-refStart:prevEnd-refStart] = saveVal
            saveStart = interval.start
            saveVal = interval.val
        prevEnd = interval.end
        interval = interval.next

    if not firstTime:
       out[saveStart-refStart:prevEnd-refStart] = saveVal

    if refStart < start:
        out[:(start-refStart)] = oob
    if end >= chromSize:
        out[(chromSize - refStart):] = oob

    lmCleanup(&lm)


def fetch(str inFile, str chrom, int start, int end, double missing=0.0, double oob=np.nan):
    """
    Fetch a BigWig interval at base pair resolution as a numpy array.

    Parameters
    ----------
    inFile : str
        Path to BigWig file.
    chrom : str
        Chromosome name.
    start : int
        Start coordinate. If start is less than zero, the beginning of the track
        is not truncated but treated as out of bounds.
    end : int
        End coordinate. If end is less than zero, the end is set to the
        chromosome size. If end is greater than the chromosome size, the end of
        the track is not truncated but treated as out of bounds.
    missing : float
        Fill-in value for unreported data in valid regions. Default is 0.
    oob : float
        Fill-in value for out-of-bounds regions. Default is NaN.

    Returns
    -------
    out : ndarray

    """
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef char *cInFile = bInFile
    cdef bytes bChrom = chrom.encode('utf-8')
    cdef char *cChrom = bChrom
    cdef bbiFile *bwf = bigWigFileOpen(cInFile)
    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chromobj
    cdef int refStart, length

    cdef np.ndarray[np.double_t, ndim=1] out

    chromList = chromobj = bbiChromList(bwf)
    while chromobj != NULL:
        if sameString(cChrom, chromobj.name):
            break
        chromobj = chromobj.next
    else:
        raise KeyError("Chromosome not found: {}".format(chrom))

    refStart = start
    if start < 0:
        start = 0
    if end < 0:
        end = chromobj.size
    length = end - refStart

    if length < 0:
        raise ValueError(
            "Interval cannot have negative length:"
            " start = {}, end = {}.".format(start, end))

    if missing == 0.0:
        out = np.zeros(length, dtype=float)
    else:
        out = np.full(length, fill_value=missing, dtype=float)

    _fetch(out, bwf, chromobj.name, start, end, refStart, chromobj.size, missing, oob)

    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bwf)

    return out


#def pileup(str infile, pd.DataFrame bed, double missing=0.0, double oob=np.nan):
#    pass
