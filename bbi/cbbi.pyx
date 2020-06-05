#!python
#cython: embedsignature=True
from __future__ import division, print_function
from six.moves.urllib.request import urlopen
from six.moves.urllib.parse import urlparse
from collections import OrderedDict
import os.path
import sys

import numpy as np
import cython

from libc.math cimport sqrt


if sys.version_info.major > 2:
    bytes_to_int = int.from_bytes
else:
    from contextlib import closing
    _urlopen = urlopen

    def urlopen(*a, **kw):
        return closing(_urlopen(*a, **kw))

    def bytes_to_int(b, byteorder):
        from struct import unpack
        fmt = '<I' if byteorder == 'little' else '>I'
        return unpack(fmt, b)[0]


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
        # if not _ucsc_may_open_url(uri):
        #     raise RuntimeError("UCSC lib cannot open this URL")
    return magic_bytes


def _check_sig(str uri):
    magic_bytes = _read_magic(uri)
    if (bytes_to_int(magic_bytes, 'little') == bigWigSig or
        bytes_to_int(magic_bytes, 'big') == bigWigSig):
        return bigWigSig
    elif (bytes_to_int(magic_bytes, 'little') == bigBedSig or
          bytes_to_int(magic_bytes, 'big') == bigBedSig):
        return bigBedSig
    else:
        return 0


def _ucsc_may_open_url(str url):
    cdef bytes bUrl = url.encode('utf-8')
    f = udcFileMayOpen(bUrl, udcDefaultDir())
    if f == NULL:
        return False
    else:
        udcFileClose(&f)
        return True


cpdef dict BBI_SUMMARY_TYPES = {
    'mean': bbiSumMean,
    'max': bbiSumMax,
    'min': bbiSumMin,
    'cov': bbiSumCoverage,
    'std': bbiSumStandardDeviation,
    'sum': bbiSumSum,
}


def is_bbi(str inFile):
    """
    Returns True if `inFile` is a path or URL to a big binary file.

    Parameters
    ----------
    inFile : str
        File path or URL

    Returns
    -------
    bool

    """
    return _check_sig(inFile) > 0


def is_bigwig(str inFile):
    """
    Returns True if `inFile` is a path or URL to a BigWig binary file.

    Parameters
    ----------
    inFile : str
        File path or URL

    Returns
    -------
    bool

    """
    return _check_sig(inFile) == bigWigSig


def is_bigbed(str inFile):
    """
    Returns True if `inFile` is a path or URL to a BigBed binary file.

    Parameters
    ----------
    inFile : str
        File path or URL

    Returns
    -------
    bool

    """
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
    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bbiFile *bbi
    if sig == bigWigSig:
        bbi = bigWigFileOpen(bInFile)
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(bInFile)
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # traverse the chromosome list
    # chromList is allocated
    cdef bbiChromInfo *chromList = bbiChromList(bbi)
    cdef bbiChromInfo *info = chromList
    cdef list c_list = []
    while info != NULL:
        c_list.append((<bytes>(info.name).decode('ascii'), info.size))
        info = info.next

    # clean up
    bbiChromInfoFreeList(&chromList)
    bbiFileClose(&bbi)

    return OrderedDict(c_list)


def zooms(str inFile):
    """
    Fetch the zoom levels of a bbi file. Returns a list of "reduction levels",
    i.e. the number of bases per summary item, i.e. the bin size.

    Parameters
    ----------
    inFile : str
        Path to BigWig or BigBed file.

    Returns
    -------
    list of int

    """
    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bbiFile *bbi
    if sig == bigWigSig:
        bbi = bigWigFileOpen(bInFile)
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(bInFile)
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # traverse the zoom list
    cdef bbiZoomLevel *zoom = bbi.levelList
    cdef list z_list = []
    while zoom != NULL:
        z_list.append(zoom.reductionLevel)
        zoom = zoom.next

    # clean up
    bbiFileClose(&bbi)

    return z_list


def info(str inFile):
    """
    Returns a dict of information about the bbi file.

    Parameters
    ----------
    inFile : str
        Path to BigWig or BigBed file.

    Returns
    -------
    dict

    """
    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bbiFile *bbi
    if sig == bigWigSig:
        bbi = bigWigFileOpen(bInFile)
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(bInFile)
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    cdef bbiSummaryElement summ = bbiTotalSummary(bbi)
    d = {
        'version': bbi.version,
        'isCompressed': bbi.uncompressBufSize > 0,
        'isSwapped': bbi.isSwapped,
        'primaryDataSize': bbi.unzoomedIndexOffset - bbi.unzoomedDataOffset,
        'zoomLevels': bbi.zoomLevels,
        'chromCount': len(chromsizes(inFile)),
        'summary': {
            'basesCovered': summ.validCount,
            'mean': summ.sumData / summ.validCount,
            'min': summ.minVal,
            'max': summ.maxVal,
            'std': sqrt(var_from_sums(summ.sumData,
                                      summ.sumSquares,
                                      summ.validCount)),
            'sum': summ.sumData,
        }
    }
    bbiFileClose(&bbi)

    return d


def fetch_intervals(
        str inFile,
        str chrom,
        int start,
        int end):
    """
    Return a generator that iterates over the records of data intervals in a
    bbi file overlapping a specified genomic query interval.

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
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bbiFile *bbi
    if sig == bigWigSig:
        bbi = bigWigFileOpen(bInFile)
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(bInFile)
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # find the chromosome
    cdef bytes chromName = chrom.encode('ascii')
    cdef int chromSize = bbiChromSize(bbi, chromName)
    if chromSize == 0:
        bbiFileClose(&bbi)
        raise KeyError("Chromosome not found: {}".format(chrom))

    # check the coordinates
    if end < 0:
        end = chromSize
    if start > chromSize:
        bbiFileClose(&bbi)
        raise ValueError(
            "Start exceeds the chromosome length, {}.".format(chromSize))
    cdef length = end - start
    if length < 0:
        bbiFileClose(&bbi)
        raise ValueError(
            "Interval cannot have negative length:"
            " start = {}, end = {}.".format(start, end))

    # clip the query range
    cdef int validStart = start, validEnd = end
    if start < 0:
        validStart = 0
    if end > chromSize:
        validEnd = chromSize

    # query
    # interval list is allocated out of lm
    cdef lm *lm = lmInit(0)
    cdef bbiInterval *bwInterval
    cdef bigBedInterval *bbInterval
    cdef tuple restTup
    cdef char *cRest = NULL
    cdef bytes bRest
    cdef str rest
    if sig == bigWigSig:
        bwInterval = bigWigIntervalQuery(bbi, chromName, validStart, validEnd, lm)
        while bwInterval != NULL:
            yield (chrom, bwInterval.start, bwInterval.end, bwInterval.val)
            bwInterval = bwInterval.next
    else:
        bbInterval = bigBedIntervalQuery(bbi, chromName, validStart, validEnd, 0, lm)
        while bbInterval != NULL:
            cRest = bbInterval.rest
            if cRest != NULL:
                bRest = cRest
                rest = bRest.decode('ascii')
                restTup = tuple(rest.split('\t'))
            else:
                restTup = ()
            yield (chrom, bbInterval.start, bbInterval.end) + restTup
            bbInterval = bbInterval.next

    # clean up
    lmCleanup(&lm)
    bbiFileClose(&bbi)


def fetch(
        str inFile,
        str chrom,
        int start,
        int end,
        int bins=-1,
        double missing=0.0,
        double oob=np.nan,
        str summary='mean'):
    """
    Read the signal data in a bbi file overlapping a genomic query interval
    into a numpy array.

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
        Start coordinate. If start is less than zero, the beginning of the
        track is not truncated but treated as out of bounds.
    end : int
        End coordinate. If end is less than zero, the end is set to the
        chromosome size. If end is greater than the chromosome size, the end of
        the track is not truncated but treated as out of bounds.
    bins : int, optional
        Number of bins to divide the query interval into for coarsegraining.
        Default (-1) means no summarization (i.e., 1 bp bins).
    missing : float, optional
        Fill-in value for unreported data in valid regions. Default is 0.
    oob : float, optional
        Fill-in value for out-of-bounds regions. Default is NaN.
    summary : str, optional
        Summary statistic to use if summarizing. Options are 'mean', 'min',
        'max', 'cov' (coverage), and 'std' (standard deviation). Default is
        'mean'.

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
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bbiFile *bbi
    cdef BbiFetchIntervals fetcher
    if sig == bigWigSig:
        bbi = bigWigFileOpen(bInFile)
        fetcher = bigWigIntervalQuery
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(bInFile)
        fetcher = bigBedCoverageIntervals
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # find the chromosome
    cdef bytes chromName = chrom.encode('ascii')
    cdef int chromSize = bbiChromSize(bbi, chromName)
    if chromSize == 0:
        bbiFileClose(&bbi)
        raise KeyError("Chromosome not found: {}".format(chrom))

    # check the coordinates
    if end < 0:
        end = chromSize
    if start > chromSize:
        bbiFileClose(&bbi)
        raise ValueError(
            "Start exceeds the chromosome length, {}.".format(chromSize))
    cdef length = end - start
    if length < 0:
        bbiFileClose(&bbi)
        raise ValueError(
            "Interval cannot have negative length:"
            " start = {}, end = {}.".format(start, end))

    # prepare the output
    cdef boolean is_summary = True
    if bins < 1:
        is_summary = False
        bins = length

    cdef np.ndarray[np.double_t, ndim=1] out = np.empty(bins, dtype=float)
    out[:] = missing

    # query
    cdef bbiSummaryType summary_type
    if is_summary:
        try:
            summary_type = BBI_SUMMARY_TYPES[summary]
        except KeyError:
            bbiFileClose(&bbi)
            raise ValueError(
                'Invalid summary type "{}". Must be one of: {}.'.format(
                    summary,
                    set(BBI_SUMMARY_TYPES.keys())))
        array_query_summarized(
            out, bins, bbi, fetcher,
            chromName, start, end, chromSize, oob, summary_type)
    else:
        array_query_full(
            out, bins, bbi, fetcher,
            chromName, start, end, chromSize, oob)

    # clean up
    bbiFileClose(&bbi)

    return out


def stackup(
        str inFile,
        chroms,
        starts,
        ends,
        int bins=-1,
        double missing=0.0,
        double oob=np.nan,
        str summary='mean'):
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
        Number of bins to summarize the data. Default (-1) means no
        aggregation.
    missing : float
        Fill-in value for unreported data in valid regions. Default is 0.
    oob : float
        Fill-in value for out-of-bounds regions. Default is NaN.
    summary : str, optional
        Summary statistic to use if summarizing. Options are 'mean', 'min',
        'max', 'cov' (coverage), and 'std' (standard deviation). Default is
        'mean'.

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
    if len(chroms_) != len(starts_) or len(starts_) != len(ends_):
        raise ValueError(
            "`chroms`, `starts`, and `ends` must have the same length")

    # open the file
    cdef bits32 sig = _check_sig(inFile)
    cdef bytes bInFile = inFile.encode('utf-8')
    cdef bbiFile *bbi
    cdef BbiFetchIntervals fetcher
    if sig == bigWigSig:
        bbi = bigWigFileOpen(bInFile)
        fetcher = bigWigIntervalQuery
    elif sig == bigBedSig:
        bbi = bigBedFileOpen(bInFile)
        fetcher = bigBedCoverageIntervals
    else:
        raise OSError("Not a bbi file: {}".format(inFile))

    # check the coordinate inputs
    if bins < 0 and len(np.unique(ends_ - starts_)) != 1:
        bbiFileClose(&bbi)
        raise ValueError(
            "Query windows must have equal size if `bins` is not specified."
        )

    # prepare output
    cdef int length = ends_[0] - starts_[0]
    cdef int n = chroms_.shape[0]
    if bins < 1:
        bins = length
    cdef np.ndarray[np.double_t, ndim=2] out = np.empty((n, bins), dtype=float)
    out[:, :] = missing

    # query
    cdef bytes chromName
    cdef int chromSize
    cdef int start, end
    cdef Py_ssize_t i
    cdef bbiSummaryType summary_type
    for i in range(n):
        chromName = chroms_[i].encode('ascii')
        chromSize = bbiChromSize(bbi, chromName)
        if chromSize == 0:
            bbiFileClose(&bbi)
            raise KeyError("Chromosome not found: {}".format(chroms_[i]))
        start = starts_[i]
        end = ends_[i]
        if start > chromSize:
            bbiFileClose(&bbi)
            raise ValueError(
                "Start exceeds the chromosome length, {}.".format(chromSize))
        if bins < 1:
            array_query_full(
                out[i, :], bins, bbi, fetcher,
                chromName, start, end, chromSize, oob)
        else:
            try:
                summary_type = BBI_SUMMARY_TYPES[summary]
            except KeyError:
                bbiFileClose(&bbi)
                raise ValueError(
                    'Invalid summary type "{}". Must be one of: {}.'.format(
                        summary,
                        set(BBI_SUMMARY_TYPES.keys())))
            array_query_summarized(
                out[i, :], bins, bbi, fetcher,
                chromName, start, end, chromSize, oob, summary_type)

    # clean up
    bbiFileClose(&bbi)

    return out


cdef inline void array_query_full(
        np.ndarray[np.double_t, ndim=1] out,
        int nbins,
        bbiFile *bbi,
        BbiFetchIntervals fetchIntervals,
        bytes chromName,
        int start,
        int end,
        int chromSize,
        double oob):

    # Clip the query range
    cdef int validStart = start, validEnd = end
    if start < 0:
        validStart = 0
    if end > chromSize:
        validEnd = chromSize

    # Fill valid regions
    # intervalList is allocated out of lm
    cdef lm *lm = lmInit(0)
    cdef bbiInterval *intervalList = fetchIntervals(
        bbi, chromName, validStart, validEnd, lm)

    cdef:
        boolean firstTime = True
        int saveStart = -1
        int prevEnd = -1
        double saveVal = -1.0
        bbiInterval *interval = intervalList
    while interval != NULL:
        if firstTime:
            saveStart = interval.start
            saveVal = interval.val
            firstTime = False
        elif not ((interval.start == prevEnd) and (interval.val == saveVal)):
            out[saveStart-start:prevEnd-start] = saveVal
            saveStart = interval.start
            saveVal = interval.val
        prevEnd = interval.end
        interval = interval.next
    if not firstTime:
        out[saveStart-start:prevEnd-start] = saveVal

    # Fill out-of-bounds regions
    if start < validStart:
        out[:(validStart - start)] = oob
    if end >= validEnd:
        out[(validEnd - start):] = oob

    lmCleanup(&lm)


cdef double var_from_sums(double sum, double sumSquares, bits64 n):
    cdef double var = sumSquares - sum*sum/n
    if n > 1:
        var /= n - 1
    return var


cdef inline void array_query_summarized(
        np.ndarray[np.double_t, ndim=1] out,
        int nbins,
        bbiFile *bbi,
        BbiFetchIntervals fetchIntervals,
        bytes chromName,
        int start,
        int end,
        int chromSize,
        double oob,
        bbiSummaryType summaryType):

    # Clip the query range
    cdef int validStart = start, validEnd = end
    if start < 0:
        validStart = 0
    if end > chromSize:
        validEnd = chromSize

    # Get the closest zoom level less than what we're looking for
    cdef int baseSize = end - start
    cdef int stepSize = baseSize // nbins
    cdef int zoomLevel = stepSize // 2
    if zoomLevel < 0:
        zoomLevel = 0
    cdef bbiZoomLevel *zoomObj = bbiBestZoom(bbi.levelList, zoomLevel)

    # Create and populate summary elements
    # elements is allocated
    cdef boolean result = False
    cdef bbiSummaryElement *elements
    AllocArray(elements, nbins)
    if zoomObj != NULL:
        result = _bbiSummariesFromZoom(
            bbi, zoomObj,
            chromName, start, end, validStart, validEnd,
            elements, nbins)
    else:
        result = _bbiSummariesFromFull(
            bbi, fetchIntervals,
            chromName, start, end, validStart, validEnd,
            elements, nbins)

    # Fill output array
    cdef double covFactor = <double>nbins / (end - start)
    cdef bbiSummaryElement *el
    cdef double val
    cdef int loc, i
    if result:
        for i in range(nbins):
            loc = start + i*stepSize
            if loc < validStart or loc >= validEnd:
                out[i] = oob
            else:
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
                        val = sqrt(var_from_sums(el.sumData,
                                                 el.sumSquares,
                                                 el.validCount))
                    elif summaryType == bbiSumSum:
                        val = el.sumData
                    else:
                        raise RuntimeError
                    out[i] = val

    # Destroy summary elements
    freeMem(elements)


cdef boolean _bbiSummariesFromZoom(
       bbiFile *bbi,
       bbiZoomLevel *zoom,
       bytes chromName,
       int start,
       int end,
       int validStart,
       int validEnd,
       bbiSummaryElement *elements,
       int nbins):
    # Look up region in index and get data at given zoom level.
    # Summarize this data in the summary array.

    cdef int chromId = bbiChromId(bbi, chromName)
    if chromId < 0:
        return False

    # Find appropriate zoom-level summary data
    # summList is allocated
    cdef bbiSummary *summList = bbiSummariesInRegion(
        zoom, bbi, chromId, validStart, validEnd)

    # Interpolate the summary data
    cdef boolean result = False
    cdef np.int64_t baseCount = end - start
    cdef np.int64_t baseStart = start
    cdef np.int64_t baseEnd
    cdef int i
    cdef bbiSummary *summ = summList
    if summ != NULL:
        for i in range(nbins):
            # Calculate end of this part of summary
            baseEnd = start + baseCount*(i+1) // nbins

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


cdef boolean _bbiSummariesFromFull(
        bbiFile *bbi,
        BbiFetchIntervals fetchIntervals,
        bytes chromName,
        int start,
        int end,
        int validStart,
        int validEnd,
        bbiSummaryElement *elements,
        int nbins):
    # Summarize data, not using zoom. Updates the summary elements.

    # Find appropriate interval elements
    # intervalList is allocated out of lm
    cdef lm *lm = lmInit(0)
    cdef bbiInterval *intervalList = NULL
    intervalList = fetchIntervals(bbi, chromName, validStart, validEnd, lm)

    # Interpolate the interval data
    cdef boolean result = False
    cdef np.int64_t baseCount = end - start
    cdef np.int64_t baseStart = start
    cdef np.int64_t baseEnd
    cdef int end1
    cdef int i
    cdef bbiInterval *interval = intervalList
    if interval != NULL:
        for i in range(nbins):
            # Calculate end of this part of summary
            baseEnd = start + baseCount*(i+1) // nbins
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
