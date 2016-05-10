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


def bigWigToBedGraph(inFile, char *chrom_, int start_, int end_):
    cdef bbiFile *bwf = bigWigFileOpen(inFile)
    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chrom
    cdef bbiInterval *intervalList
    cdef bbiInterval *interval
    cdef lm *lm
    
    cdef boolean firstTime
    cdef int saveStart, prevEnd
    cdef double saveVal

    chromList = chrom = bbiChromList(bwf)
    while chrom != NULL:
        if sameString(chrom_, chrom.name):
            # we found the right chromosome
            lm = lmInit(0)
            chromName = chrom.name
            start = 0
            end = chrom.size
            if start_ > 0:
                start = start_
            if end_ > 0:
                end = end_

            # fetch the intervals
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
                    print(chromName, saveStart, prevEnd, saveVal)
                    saveStart = interval.start
                    saveVal = interval.val
                prevEnd = interval.end
                interval = interval.next
            if not firstTime:
               print(chromName, saveStart, prevEnd, saveVal)
            lmCleanup(&lm)
            break
        else:
            chrom = chrom.next

    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bwf)


def fetch(str inFile, str chrom, int start, int end, double missing=0.0, double oob=np.nan):
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bytes bChrom = chrom.encode('utf-8')
    cdef char *cInFile = bInFile
    cdef char *cChrom = bChrom
    cdef bbiFile *bwf = bigWigFileOpen(cInFile)

    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chromobj
    cdef bbiInterval *intervalList
    cdef bbiInterval *interval
    cdef lm *lm
    cdef boolean firstTime
    cdef int saveStart, prevEnd, refStart, length
    cdef double saveVal
    cdef np.ndarray[np.double_t, ndim=1] out

    chromList = chromobj = bbiChromList(bwf)
    while chromobj != NULL:
        if sameString(cChrom, chromobj.name):
            # we found the right chromosome
            lm = lmInit(0)
            chromName = chromobj.name
            refStart = start
            if start <= 0:
                start = 0
            if end <= 0:
                end = chromobj.size

            length = end - refStart
            if missing == 0:
                out = np.zeros(length, dtype=float)
            else:
                out = np.full(length, fill_value=missing, dtype=float)

            # fetch the intervals
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
            if end >= chromobj.size:
                out[(chromobj.size - refStart):] = oob

            lmCleanup(&lm)
            break
        else:
            chromobj = chromobj.next

    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bwf)

    return out
