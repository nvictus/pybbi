#!python
#cython: embedsignature=True
from __future__ import division, print_function
import numpy as np
import cython


def fetch(str inFile,
        str chrom, 
        int start, 
        int end, 
        int bins=-1, 
        double missing=0.0, 
        double oob=np.nan):
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
    bins : int
        Number of bins to summarize the data. Default (-1) means no aggregation.
    missing : float
        Fill-in value for unreported data in valid regions. Default is 0.
    oob : float
        Fill-in value for out-of-bounds regions. Default is NaN.

    Returns
    -------
    out : 1D ndarray

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
    
    if bins < 1:
        bins = length
    
    if missing == 0.0:
        out = np.zeros(bins, dtype=float)
    else:
        out = np.full(bins, fill_value=missing, dtype=float)

    if bins < 1:
        bwFetchFull(out, bwf, chromobj.name, start, end, refStart,
            chromobj.size, missing, oob)
    else:
        bwFetchSummarized(out, bwf, chromobj.name, start, end, refStart,
            chromobj.size, bins, missing, oob)

    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bwf)

    return out


def stackup(str inFile, 
        chroms,
        starts,
        ends,
        int bins=-1,
        double missing=0.0,
        double oob=np.nan):
    """
    Vertically stack equal-length bigwig intervals.

    Parameters
    ----------
    inFile : str
        Path to BigWig file.
    chrom : array-like of str
        Chromosome names.
    start : array-like of int
        Start coordinates. If start is less than zero, the beginning of the
        track is not truncated but treated as out of bounds.
    end : array-like of int
        End coordinates. If end is less than zero, the end is set to the
        chromosome size. If end is greater than the chromosome size, the end of
        the track is not truncated but treated as out of bounds.
    bins : int
        Number of bins to summarize the data. Default (-1) means no aggregation.
    missing : float
        Fill-in value for unreported data in valid regions. Default is 0.
    oob : float
        Fill-in value for out-of-bounds regions. Default is NaN.

    Returns
    -------
    out : ndarray

    """
    cdef np.ndarray[object, ndim=1] chroms_ = np.asarray(chroms, dtype=object)
    cdef np.ndarray[np.int_t, ndim=1] starts_ = np.asarray(starts, dtype=int)
    cdef np.ndarray[np.int_t, ndim=1] ends_ = np.asarray(ends, dtype=int)
    
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef char *cInFile = bInFile
    cdef bbiFile *bwf = bigWigFileOpen(cInFile)
    
    cdef int n = chroms_.shape[0]
    cdef length = ends_[0] - starts_[0]

    cdef int start, end, refStart, i
    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chromobj
    cdef bytes bChrom
    cdef char *cChrom

    if not len(np.unique(ends_ - starts_)) == 1:
        raise ValueError("Windows must have equal size")

    if bins < 1:
        bins = length
        
    if missing == 0.0:
        out = np.zeros((n, bins), dtype=float)
    else:
        out = np.full((n, bins), fill_value=missing, dtype=float)

    for i in range(n):
        bChrom = chroms_[i].encode('utf-8')
        cChrom = bChrom
        chromList = chromobj = bbiChromList(bwf)
        while chromobj != NULL:
            if sameString(cChrom, chromobj.name):
                break
            chromobj = chromobj.next
        else:
            raise KeyError(
                "Chromosome not found: {}".format(bChrom.decode('utf-8')))

        start = starts_[i]
        end = ends_[i]
        refStart = start
        if start < 0:
            start = 0

        if bins < 1:
            bwFetchFull(out[i, :], bwf, chromobj.name, start, end, refStart,
                chromobj.size, missing, oob)
        else:
            bwFetchSummarized(out[i, :], bwf, chromobj.name, start, end,
                refStart, chromobj.size, bins, missing, oob)
            
    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bwf)

    return out


cdef inline void bwFetchFull(
        np.ndarray[np.double_t, ndim=1] out, 
        bbiFile *bwf, 
        char *chromName, 
        int start,
        int end,
        int refStart,
        int chromSize,
        double missing, 
        double oob):
    
    cdef bbiInterval *intervalList
    cdef bbiInterval *interval
    cdef lm *lm
    cdef boolean firstTime
    cdef int saveStart, prevEnd
    cdef double saveVal

    lm = lmInit(0)
    intervalList = interval = bigWigIntervalQuery(
        bwf, chromName, start, end, lm)
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


cdef inline void bwFetchSummarized(
        np.ndarray[np.double_t, ndim=1] out, 
        bbiFile *bwf, 
        char *chromName, 
        int start,
        int end,
        int refStart,
        int chromSize,
        int nbins,
        double missing, 
        double oob):
    
    # Get the closest zoom level less than what we're looking for
    cdef int baseSize = end - refStart
    cdef int fullReduction = baseSize // nbins
    cdef int zoomLevel = fullReduction // 2
    if zoomLevel < 0:
        zoomLevel = 0
    cdef bbiZoomLevel *zoomObj = bbiBestZoom(bwf.levelList, zoomLevel)

    # Create summary elements
    cdef bbiSummaryElement *elements
    AllocArray(elements, nbins)

    # Populate summary elements
    cdef boolean result = False
    if zoomObj != NULL:
        result = _bbiSummarizeFromZoom(zoomObj, bwf, chromName, start, end,
            refStart, nbins, elements)
    else:
        result = _bbiSummarizeFromFull(bwf, chromName, start, end, refStart,
            nbins, bigWigIntervalQuery, elements)

    # Fill array
    cdef double covFactor = <double>nbins / (end - start)
    cdef bbiSummaryElement *el
    cdef double val
    cdef summaryType = bbiSumMean
    cdef int i
    if result:
        for i in range(nbins):
            el = &elements[i]
            if el.validCount > 0:
                if summaryType == bbiSumMean:
                    val = el.sumData / el.validCount
                elif summaryType == bbiSumMax:
                    val = el.maxVal
                elif summaryType == bbiSumMin:
                    val = el.minVal
                elif summaryType == bbiSumCoverage:
                    val = covFactor * el.validCount
                elif summaryType == bbiSumStandardDeviation:
                    val = 0.0 #calcStdFromSums(el.sumData, el.sumSquares, el.validCount)
                else:
                    raise RuntimeError
            out[i] = val

    # Destroy summary element array
    freeMem(elements)    


cdef boolean _bbiSummarizeFromZoom(
       bbiZoomLevel *zoom,
       bbiFile *bbi, 
       char *chrom, 
       int start, 
       int end,
       int refStart,
       int nbins, 
       bbiSummaryElement *elements):
    # Look up region in index and get data at given zoom level.
    # Summarize this data in the summary array.

    cdef int chromId = bbiChromId(bbi, chrom)
    if chromId < 0:
        return False

    # Find appropriate zoom-level summary data
    cdef bbiSummary *summList = bbiSummariesInRegion(
        zoom, bbi, chromId, start, end)
    
    # Interpolate the summary data
    cdef boolean result = False
    cdef np.int64_t baseCount = end - refStart
    cdef np.int64_t baseStart = refStart
    cdef np.int64_t baseEnd
    cdef int i
    cdef bbiSummary *summ = summList
    if summ != NULL:
        for i in range(nbins):
            # Calculate end of this part of summary
            baseEnd = refStart + baseCount*(i+1) // nbins

            # Advance summ to skip over parts we are no longer interested in.
            while summ != NULL and summ.end <= baseStart:
                summ = summ.next

            # Update element[i]
            if bbiSummarySlice(bbi, baseStart, baseEnd, summ, &elements[i]):
                result = True

            # Next time round start where we left off.
            baseStart = baseEnd

        slFreeList(&summList)

    return result


cdef boolean _bbiSummarizeFromFull(
        bbiFile *bbi, 
        char *chrom, 
        int start, 
        int end,
        int refStart,
        int nbins,
        BbiFetchIntervals fetchIntervals,
        bbiSummaryElement *elements):
    # Summarize data, not using zoom. Updates the summary elements.

    # Find appropriate interval elements
    cdef lm *lm = lmInit(0)
    cdef bbiInterval *intervalList = NULL
    intervalList = fetchIntervals(bbi, chrom, start, end, lm)
    
    # Interpolate the interval data
    cdef boolean result = False
    cdef np.int64_t baseCount = end - refStart
    cdef np.int64_t baseStart = refStart
    cdef np.int64_t baseEnd
    cdef int end1
    cdef int i
    cdef bbiInterval *interval = intervalList
    if interval != NULL:
        for i in range(nbins):
            # Calculate end of this part of summary
            baseEnd = refStart + baseCount*(i+1) // nbins
            end1 = baseEnd
            if end1 == baseStart:
                end1 = baseStart + 1

            # Advance interval to skip over parts we are no longer interested in.
            while interval != NULL  and interval.end <= baseStart:
                interval = interval.next

            # Update element[i]
            if bbiIntervalSlice(bbi, baseStart, end1, interval, &elements[i]):
                result = True

            # Next time round start where we left off.
            baseStart = baseEnd
        
    lmCleanup(&lm)
    
    return result
