#!python
#cython: embedsignature=True
from __future__ import division, print_function
from six.moves.urllib.request import urlopen
from six.moves.urllib.parse import urlparse
from collections import OrderedDict
import os.path

import numpy as np
import cython


cdef inline char *cstring_encode(str text, str encoding='ascii'):
    cdef bytes bText = text.encode(encoding)
    cdef char *cText = bText
    return cText


cdef inline str cstring_decode(char *cText, str encoding='ascii'):
    cdef bytes bText = cText
    cdef str text = bText.decode(encoding)
    return text


def _is_url(str uri):
    return urlparse(uri).scheme != ""


def _read_magic(str uri):
    cdef bytes magic_bytes
    if not _is_url(uri):
        if not os.path.isfile(uri):
            raise OSError("File not found: {}".format(uri))
        with open(uri, 'rb') as f:
            magic_bytes = f.read(4)
    else:
        with urlopen(uri) as r:
            code = r.getcode()
            if code >= 400:
                raise OSError("Status {}: Couldn't open {}".format(code, uri))
            magic_bytes = r.read(4)
    return magic_bytes


def _check_sig(str uri):
    magic_bytes = _read_magic(uri)
    if (int.from_bytes(magic_bytes, 'little') == bigWigSig or
        int.from_bytes(magic_bytes, 'big') == bigWigSig):
        return bigWigSig
    elif (int.from_bytes(magic_bytes, 'little') == bigBedSig or
          int.from_bytes(magic_bytes, 'big') == bigBedSig):
        return bigBedSig
    else:
        return 0


def is_bbi(str inFile):
    return _check_sig(inFile) > 0


def is_bigwig(str inFile):
    return _check_sig(inFile) == bigWigSig


def is_bigbed(str inFile):
    return _check_sig(inFile) == bigBedSig


def chromsizes(str inFile):
    """
    Fetch the chromosome list of a bbi file. Returns an ordered dictionary of
    chromosome names mapped to their sizes in bp.

    Parameters
    ----------
    inFile : str
        Path to BigWig or BigBed file.

    Returns
    -------
    OrderedDict (str -> int)

    """
    cdef bits32 sig = _check_sig(inFile)
    cdef char *cInFile = cstring_encode(inFile, 'utf-8')
    cdef bbiFile *bbi
    cdef BbiFetchIntervals fetcher
    if sig == bigWigSig:
        bbi = bigWigFileOpen(cInFile)
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(cInFile)
    else:
        raise OSError("Not a bbi file: {}".format(inFile))
    cdef bbiChromInfo *chromList = bbiChromList(bbi)
    cdef bbiChromInfo *info = chromList 
    cdef list c_list = []
    while info != NULL:
        c_list.append( (cstring_decode(info.name, 'ascii'), info.size) )
        info = info.next
    bbiChromInfoFreeList(&chromList)
    return OrderedDict(c_list)


def fetch_intervals(
        str inFile,
        str chrom, 
        int start, 
        int end):
    """
    Return a generator that iterates over the records of data intervals in a bbi
    file overlapping a specified genomic query interval.

    Parameters
    ----------
    inFile : str
        Path to BigWig or BigBed file.
    chrom : str
        Chromosome name.
    start : int
        Start coordinate.
    end : int
        End coordinate. If end is less than zero, the end is set to the
        chromosome size.

    Yields
    ------
    tuple
        bedGraph or BED record

    """
    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef char *cInFile = cstring_encode(inFile, 'utf-8')
    cdef bbiFile *bbi
    cdef BbiFetchIntervals fetcher
    if sig == bigWigSig:
        bbi = bigWigFileOpen(cInFile)
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(cInFile)
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # find the chromosome
    cdef char *cChrom = cstring_encode(chrom, 'ascii')
    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chromobj
    chromList = chromobj = bbiChromList(bbi)
    while chromobj != NULL:
        if sameString(cChrom, chromobj.name):
            break
        chromobj = chromobj.next
    else:
        raise KeyError("Chromosome not found: {}".format(chrom))

    # check the coordinate inputs
    cdef int refStart, length
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

    # query
    cdef lm *lm = lmInit(0)
    cdef bbiInterval *bwInterval
    cdef bigBedInterval *bbInterval
    cdef tuple restTup
    cdef char *cRest
    cdef str rest
    if sig == bigWigSig:
        bwInterval = bigWigIntervalQuery(bbi, cChrom, start, end, lm)
        while bwInterval != NULL:
            yield (chrom, bwInterval.start, bwInterval.end, bwInterval.val)
            bwInterval = bwInterval.next
    else:
        bbInterval = bigBedIntervalQuery(bbi, cChrom, start, end, 0, lm)
        while bbInterval != NULL:
            cRest = bbInterval.rest
            if cRest != NULL:
                rest = cstring_decode(cRest, 'ascii')
                restTup = tuple(rest.split('\t'))
            else:
                restTup = ()
            yield (chrom, bbInterval.start, bbInterval.end) + restTup
            bbInterval = bbInterval.next

    # clean up
    lmCleanup(&lm)
    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bbi)


def fetch(
        str inFile,
        str chrom, 
        int start, 
        int end, 
        int bins=-1, 
        double missing=0.0, 
        double oob=np.nan):
    """
    Read the signal data in a bbi file overlapping a genomic query interval into
    a numpy array. 

    If a number of bins is requested, this will interpolate the file's stored
    "summary" data at the closest available zoom level. Otherwise, the data
    is returned at base pair resolution (default).

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
        Number of bins to divide the query interval into for coarse graining. 
        Default (-1) means no coarse graining (i.e., 1 bp bins).
    missing : float
        Fill-in value for unreported data in valid regions. Default is 0.
    oob : float
        Fill-in value for out-of-bounds regions. Default is NaN.

    Returns
    -------
    out : 1D ndarray

    Notes
    -----
    * The signal data for a BigBed file is its coverage.
    * Currently, the only summary statistic supported is 'mean'.

    """
    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef char *cInFile = cstring_encode(inFile, 'utf-8')
    cdef bbiFile *bbi
    cdef BbiFetchIntervals fetcher
    if sig == bigWigSig:
        bbi = bigWigFileOpen(cInFile)
        fetcher = bigWigIntervalQuery
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(cInFile)
        fetcher = bigBedCoverageIntervals
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # find the chromosome
    cdef char *cChrom = cstring_encode(chrom, 'ascii')
    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chromobj
    chromList = chromobj = bbiChromList(bbi)
    while chromobj != NULL:
        if sameString(cChrom, chromobj.name):
            break
        chromobj = chromobj.next
    else:
        raise KeyError("Chromosome not found: {}".format(chrom))

    # check the coordinate inputs
    cdef int refStart, length
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
    
    # prepare the output
    cdef np.ndarray[np.double_t, ndim=1] out

    if bins < 1:
        bins = length

    if missing == 0.0:
        out = np.zeros(bins, dtype=float)
    else:
        out = np.full(bins, fill_value=missing, dtype=float)

    if bins < 1:
        _bbiFetchFull(
            out, bbi, fetcher, 
            chromobj.name, start, end, refStart, chromobj.size, 
            missing, oob)
    else:
        _bbiFetchSummarized(
            out, bbi, fetcher, 
            chromobj.name, start, end, refStart, chromobj.size, 
            bins, missing, oob)

    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bbi)

    return out


def stackup(
        str inFile, 
        chroms,
        starts,
        ends,
        int bins=-1,
        double missing=0.0,
        double oob=np.nan):
    """
    Vertically stack signal tracks from equal-length bbi query intervals.

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
    out : 2D ndarray

    See Also
    --------
    fetch : Fetch the signal track of a single interval

    """
    cdef np.ndarray[object, ndim=1] chroms_ = np.asarray(chroms, dtype=object)
    cdef np.ndarray[np.int_t, ndim=1] starts_ = np.asarray(starts, dtype=int)
    cdef np.ndarray[np.int_t, ndim=1] ends_ = np.asarray(ends, dtype=int)
    
    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef char *cInFile = cstring_encode(inFile, 'utf-8')
    cdef bbiFile *bbi
    cdef BbiFetchIntervals fetcher
    if sig == bigWigSig:
        bbi = bigWigFileOpen(cInFile)
        fetcher = bigWigIntervalQuery
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(cInFile)
        fetcher = bigBedCoverageIntervals
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # check the coordinate inputs
    if not len(np.unique(ends_ - starts_)) == 1:
        raise ValueError("Windows must have equal size")
    
    # prepare output
    cdef int length = ends_[0] - starts_[0]
    cdef int n = chroms_.shape[0]
    if bins < 1:
        bins = length
    if missing == 0.0:
        out = np.zeros((n, bins), dtype=float)
    else:
        out = np.full((n, bins), fill_value=missing, dtype=float)

    # query
    cdef int start, end, refStart, i
    cdef bbiChromInfo *chromList
    cdef bbiChromInfo *chromobj
    cdef char *cChrom
    for i in range(n):
        cChrom = cstring_encode(chroms_[i], 'ascii')
        chromList = chromobj = bbiChromList(bbi)
        while chromobj != NULL:
            if sameString(cChrom, chromobj.name):
                break
            chromobj = chromobj.next
        else:
            raise KeyError(
                "Chromosome not found: {}".format(chroms_[i]))

        start = starts_[i]
        end = ends_[i]
        refStart = start
        if start < 0:
            start = 0

        if bins < 1:
            _bbiFetchFull(
                out[i, :], bbi, fetcher, 
                chromobj.name, start, end, refStart, chromobj.size, 
                missing, oob)
        else:
            _bbiFetchSummarized(
                out[i, :], bbi, fetcher, 
                chromobj.name, start, end, refStart, chromobj.size, 
                bins, missing, oob)
            
    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bbi)

    return out


cdef inline void _bbiFetchFull(
        np.ndarray[np.double_t, ndim=1] out, 
        bbiFile *bwf, 
        BbiFetchIntervals fetchIntervals,
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
    intervalList = interval = fetchIntervals(bwf, chromName, start, end, lm)
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


cdef inline void _bbiFetchSummarized(
        np.ndarray[np.double_t, ndim=1] out, 
        bbiFile *bwf, 
        BbiFetchIntervals fetchIntervals,
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
        result = _bbiSummarizedFromZoom(zoomObj, bwf, chromName, start, end,
            refStart, nbins, elements)
    else:
        result = _bbiSummarizedFromFull(bwf, chromName, start, end, refStart,
            nbins, fetchIntervals, elements)

    # Fill output array
    cdef summaryType = bbiSumMean  # Just support mean for now
    cdef double covFactor = <double>nbins / (end - start)
    cdef bbiSummaryElement *el
    cdef double val
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


cdef boolean _bbiSummarizedFromZoom(
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


cdef boolean _bbiSummarizedFromFull(
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
