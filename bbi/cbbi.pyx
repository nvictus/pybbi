#!python
#cython: embedsignature=True
from urllib.request import urlopen
from urllib.parse import urlparse
import io
import os.path as op

import numpy as np

from libc.math cimport sqrt
from .cbbi cimport (
    asObject, bbiSumMean, bbiSumMax, bbiSumMin, bbiSumCoverage,
    bbiSumStandardDeviation, bbiSumSum
)

bytes_to_int = int.from_bytes

cdef dict BBI_SUMMARY_TYPES = {
    'mean': bbiSumMean,
    'max': bbiSumMax,
    'min': bbiSumMin,
    'cov': bbiSumCoverage,
    'std': bbiSumStandardDeviation,
    'sum': bbiSumSum,
}

# Map AutoSql types to pandas-compatible dtypes
# http://genomewiki.ucsc.edu/index.php/AutoSql
cdef dict AUTOSQL_TYPE_MAP = {
    'string': 'object',   # varchar(255)
    'char': 'object',     # char
    'double': 'float64',  # double precision floating point.
    'float': 'float32',   # single precision floating point.
    'int': 'int32',       # signed 32 bit integer
    'uint': 'uint32',     # unsigned 32 bit integer
    'short': 'int16',     # signed 16 bit integer
    'ushort': 'uint16',   # unsigned 16 bit integer
    'byte': 'int8',       # signed 8 bit integer
    'ubyte': 'uint8',     # unsigned 8 bit integer
    'off': 'int64',       # 64 bit integer
    'bigint': 'int64',    # 64 bit integer

    # Exotic types will be treated as object for now
    'lstring': 'object',  # longblob
    'object': 'object',   # longblob
    'simple': 'object',   # longblob
    'enum': 'object',     # enum?
    'set': 'object',      # set?
}


def _is_url(str uri):
    return urlparse(uri).scheme != ""


def _ucsc_may_open_url(str url):
    cdef bytes bUrl = url.encode('utf-8')
    f = udcFileMayOpen(bUrl, udcDefaultDir())
    if f == NULL:
        return False
    else:
        udcFileClose(&f)
        return True


def _read_magic(str uri):
    cdef bytes magic_bytes
    if not _is_url(uri):
        if not op.isfile(uri):
            raise OSError("File not found: {}".format(uri))
        with io.open(uri, 'rb') as f:
            magic_bytes = f.read(4)
    else:
        with urlopen(uri) as r:
            code = r.getcode()
            if code >= 400:
                raise OSError("Status {}: Couldn't open {}".format(code, uri))
            magic_bytes = r.read(4)
    return magic_bytes


def _read_sig(str uri):
    magic_bytes = _read_magic(uri)
    if (bytes_to_int(magic_bytes, 'little') == bigWigSig or
        bytes_to_int(magic_bytes, 'big') == bigWigSig):
        return bigWigSig
    elif (bytes_to_int(magic_bytes, 'little') == bigBedSig or
          bytes_to_int(magic_bytes, 'big') == bigBedSig):
        return bigBedSig
    else:
        return 0


def _check_sig(str uri):
    cdef bits32 sig = _read_sig(uri)
    if sig != bigWigSig and sig != bigBedSig:
        raise OSError("Not a bbi file: {}".format(uri))
    return sig


cdef double var_from_sums(double sum, double sumSquares, bits64 n):
    cdef double var = sumSquares - sum*sum/n
    if n > 1:
        var /= n - 1
    return var


def is_bbi(inFile):
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
    return _read_sig(inFile) > 0


def is_bigwig(inFile):
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
    return _read_sig(inFile) == bigWigSig


def is_bigbed(inFile):
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
    return _read_sig(inFile) == bigBedSig


def open(str inFile):
    """
    Open a big binary file.

    Parameters
    ----------
    inFile : str
        File path or URL

    Returns
    -------
    BBIFile
    """
    return BBIFile(inFile)


cdef class BBIFile:
    """
    Interface to a UCSC Big Binary Indexed (BBI) file.

    The resource may be a bigWig or a bigBed file.
    BigBed AutoSql schemas are supported.
    """
    cdef bbiFile *bbi
    cdef bits32 sig
    cdef readonly str path
    cdef readonly bint is_remote
    cdef readonly bint is_bigwig
    cdef readonly bint is_bigbed

    def __cinit__(self, str inFile):
        cdef bytes bInFile = inFile.encode('utf-8')
        self.sig = _check_sig(inFile)
        if self.sig == bigWigSig:
            self.bbi = bigWigFileOpen(bInFile)
        else:
            self.bbi = bigBedFileOpen(bInFile)
        self.path = inFile
        self.is_remote = _is_url(inFile)
        self.is_bigwig = self.sig == bigWigSig
        self.is_bigbed = self.sig == bigBedSig

    def __dealloc__(self):
        if self.bbi != NULL:
            bbiFileClose(&self.bbi)

    def __enter__(self):
        return self

    def __exit__(self, typ, val, tb):
        self.close()

    def close(self):
        if self.bbi != NULL:
            bbiFileClose(&self.bbi)

    @property
    def closed(self):
        return self.bbi == NULL

    def read_autosql(self):
        if self.bbi == NULL:
            raise OSError("File closed")

        if self.is_bigwig:
            return None

        # Try to read autosql definition string
        # Otherwise use default autosql BED definition based on number of fields
        cdef char *cText = bigBedAutoSqlText(self.bbi)
        if cText == NULL:
            cText = bedAsDef(self.bbi.definedFieldCount, self.bbi.fieldCount)
        cdef str text = (<bytes>cText).decode('ascii')
        freeMem(cText)
        return text

    @property
    def schema(self):
        if self.bbi == NULL:
            raise OSError("File closed")

        if self.is_bigwig:
            return {
                'name': 'bedGraph',
                'comment': '',
                'columns': ['chrom', 'start', 'end', 'value'],
                'dtypes': {
                    'chrom': 'object',
                    'start': 'uint32',
                    'end': 'uint32',
                    'value': 'float32',
                },
                'description': {
                    'chrom': 'Reference sequence chromosome or scaffold',
                    'start': 'Start position in chromosome',
                    'end': 'End position in chromosome',
                    'value': 'Data value',
                }
            }

        # Try to read autosql definition string
        # Otherwise use default autosql BED definition based on number of fields
        cdef char *cText = bigBedAutoSqlText(self.bbi)
        if cText == NULL:
            cText = bedAsDef(self.bbi.definedFieldCount, self.bbi.fieldCount)

        # Parse definition string into an object and free the string
        # cdef str raw_text = (<bytes>cText).decode('ascii')
        cdef asObject *o = asParseText(cText)

        # Extract metadata
        cdef str as_name = (<bytes>(o.name)).decode('ascii')
        cdef str comment = (<bytes>(o.comment)).decode('ascii')

        # Extract column definitions
        cdef asColumn *col = o.columnList
        cdef asTypeInfo *typ
        cdef str name, t_name, d_name
        cdef list columns = []
        cdef dict dtypes = {}
        cdef dict describe = {}
        while col != NULL:

            # field name
            typ = col.lowType
            name = (<bytes>(col.name)).decode('ascii')
            if name in ('chromStart', 'chromEnd'):
                name = name[5:].lower()
            columns.append(name)

            # field dtype
            if col.isList or col.isArray:
                dtypes[name] = 'object'
            elif name != 'reserved':
                t_name = (<bytes>(typ.name)).decode('ascii')
                if t_name in AUTOSQL_TYPE_MAP:
                    dtypes[name] = AUTOSQL_TYPE_MAP[t_name]
                else:
                    dtypes[name] = 'object'

            # field description
            d_name = (<bytes>(col.comment)).decode('ascii')
            describe[name] = d_name

            col = col.next

        # Clean up
        freeMem(cText)
        asObjectFree(&o)

        return {
            'name': as_name,
            'comment': comment,
            'columns': columns,
            'dtypes': dtypes,
            'description': describe,
        }

    @property
    def chromsizes(self):
        """
        An ordered dictionary of chromosome names mapped to their sizes in bp.

        """
        if self.bbi == NULL:
            raise OSError("File closed")

        # traverse the chromosome list
        cdef bbiChromInfo *chromList = bbiChromList(self.bbi)
        cdef bbiChromInfo *info = chromList
        cdef list c_list = []
        while info != NULL:
            c_list.append((<bytes>(info.name).decode('ascii'), info.size))
            info = info.next

        # clean up
        bbiChromInfoFreeList(&chromList)

        return dict(c_list)

    @property
    def zooms(self):
        """
        A list of "reduction levels" (bin sizes), i.e. the number of bases per
        summary item.

        """
        if self.bbi == NULL:
            raise OSError("File closed")

        # traverse the zoom list
        cdef bbiZoomLevel *zoom = self.bbi.levelList
        cdef list z_list = []
        while zoom != NULL:
            z_list.append(zoom.reductionLevel)
            zoom = zoom.next

        return z_list

    @property
    def info(self):
        """
        A dict of information about the bbi file.

        """
        if self.bbi == NULL:
            raise OSError("File closed")

        cdef bbiSummaryElement summ = bbiTotalSummary(self.bbi)

        return {
            'version': self.bbi.version,
            'isCompressed': self.bbi.uncompressBufSize > 0,
            'isSwapped': self.bbi.isSwapped,
            'primaryDataSize': self.bbi.unzoomedIndexOffset - self.bbi.unzoomedDataOffset,
            'zoomLevels': self.bbi.zoomLevels,
            'chromCount': len(self.chromsizes),
            'summary': {
                'basesCovered': summ.validCount,
                'sum': summ.sumData,
                'mean': summ.sumData / summ.validCount,
                'min': summ.minVal,
                'max': summ.maxVal,
                'std': sqrt(var_from_sums(summ.sumData,
                                          summ.sumSquares,
                                          summ.validCount)),
            }
        }

    def fetch(
        self,
        str chrom,
        int start,
        int end,
        int bins=-1,
        double missing=0.0,
        double oob=np.nan,
        str summary='mean',
        bint exact=False,
    ):
        """
        Read the signal data in a BBI file overlapping a genomic query interval
        into a numpy array.

        If a number of bins is requested, this will return a summary statistic
        for each bin. Otherwise, the data is returne unsummarized at base pair
        resolution (default).

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            Start coordinate. If start is less than zero, the beginning of the
            track is not truncated but treated as out of bounds.
        end : int
            End coordinate. If end is less than zero, the end is set to the
            chromosome size. If end is greater than the chromosome size, the
            end of the track is not truncated but treated as out of bounds.
        bins : int, optional
            Number of bins to divide the query interval into for coarsegraining.
            Default (-1) means no summarization (i.e., 1 bp bins).
        missing : float, optional [default: 0.0]
            Fill-in value for unreported data in valid regions. Default is 0.
        oob : float, optional [default: NaN]
            Fill-in value for out-of-bounds regions. Default is NaN.
        summary : str, optional [default: 'mean']
            Summary statistic to use if ``bins`` is specified. Options are
            'mean', 'sum', 'min', 'max', 'cov' (coverage), and 'std'
            (standard deviation). Default is 'mean'.
        exact : bool, optional [default: False]
            If True and ``bins`` is specified, return exact summary statistic
            values instead of interpolating from the optimal zoom level.
            Default is False.

        Returns
        -------
        1D ndarray

        Notes
        -----
        Values returned for a BigWig file are derived from its quantitative
        "value" field. A BigWig file encodes a step function, and the value at
        a base is given by the signal value of the unique interval that
        contains that base.

        Values returned for a BigBed file are derived from the coverage
        (i.e., pileup) of the intervals. A BigBed file encodes a collection
        of (possibly overlapping) intervals which may or may not be associated
        with quantitative scores. The "value" at given base used here
        summarizes the number of intervals overlapping that base, not a score.

        If a number of bins is requested and ``exact`` is False, the summary
        data is interpolated from the closest available zoom level. If you
        need accurate summary data and are okay with small trade-off in speed,
        set ``exact`` to True.

        See Also
        --------
        stackup : Stack signal tracks from many query intervals into a matrix
        """
        if self.bbi == NULL:
            raise OSError("File closed")

        cdef BbiFetchIntervals fetcher
        if self.is_bigwig:
            fetcher = bigWigIntervalQuery
        elif self.is_bigbed:
            fetcher = bigBedCoverageIntervals

        # find the chromosome
        cdef bytes chromName = chrom.encode('ascii')
        cdef int chromSize = bbiChromSize(self.bbi, chromName)
        if chromSize == 0:
            raise KeyError("Chromosome not found: {}".format(chrom))

        # check the coordinates
        if end < 0:
            end = chromSize
        if start > chromSize:
            raise ValueError(
                "Start exceeds the chromosome length, {}.".format(chromSize))
        cdef int length = end - start
        if length < 0:
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
                raise ValueError(
                    'Invalid summary type "{}". Must be one of: {}.'.format(
                        summary,
                        set(BBI_SUMMARY_TYPES.keys())))
            array_query_summarized(
                out, bins, self.bbi, fetcher,
                chromName, start, end, chromSize, oob, summary_type, exact)
        else:
            array_query_full(
                out, bins, self.bbi, fetcher,
                chromName, start, end, chromSize, oob)

        return out

    def stackup(
        self,
        chroms,
        starts,
        ends,
        int bins=-1,
        double missing=0.0,
        double oob=np.nan,
        str summary='mean',
        bint exact=False,
    ):
        """
        Vertically stack signal tracks from many query intervals into a matrix.

        Parameters
        ----------
        chroms : array-like of str
            Chromosome names of the query intervals.
        starts : array-like of int
            Start coordinates. If start is less than zero, the beginning of the
            track is not truncated but treated as out of bounds.
        ends : array-like of int
            End coordinates. If end is less than zero, the end is set to the
            chromosome size. If end is greater than the chromosome size, the
            end of the track is not truncated but treated as out of bounds.
        bins : int
            Number of bins to summarize the data. Default (-1) means no
            aggregation: each bin is 1 base pair and the query intervals must
            have the same length. If bins > 0, the input query intervals can
            have different lengths, and the data will be rescaled to the same
            number of bins.
        missing : float
            Fill-in value for unreported data in valid regions. Default is 0.
        oob : float
            Fill-in value for out-of-bounds regions. Default is NaN.
        summary : str, optional
            Summary statistic to use if summarizing. Options are 'mean', 'min',
            'max', 'cov' (coverage), and 'std' (standard deviation). Default is
            'mean'.
        exact : bool, optional [default: False]
            If True and ``bins`` is specified, return exact summary statistic
            values instead of interpolating from the optimal zoom level.
            Default is False.

        Returns
        -------
        out : 2D ndarray of float64 (n_intervals, n_bins)

        Notes
        -----
        If not summarizing, the input intervals must have the same length.

        If a number of bins is requested, intervals may have variable lengths
        and the data from each interval will be rescaled to the same number
        of bins (e.g., aligning gene bodies).

        If a number of bins is requested and ``exact`` is False, the summary
        data is interpolated from the closest available zoom level. If you
        need accurate summary data and are okay with small trade-off in
        speed, set ``exact`` to True.

        See Also
        --------
        fetch : Fetch the signal track of a single interval
        """
        if self.bbi == NULL:
            raise OSError("File closed")

        cdef np.ndarray[object, ndim=1] chroms_ = np.asarray(chroms, dtype=object)
        cdef np.ndarray[np.int64_t, ndim=1] starts_ = np.asarray(starts, dtype=int)
        cdef np.ndarray[np.int64_t, ndim=1] ends_ = np.asarray(ends, dtype=int)
        if bins < 0 and (len(chroms_) != len(starts_) or len(starts_) != len(ends_)):
            raise ValueError(
                "`chroms`, `starts`, and `ends` must have the same length"
            )

        cdef BbiFetchIntervals fetcher
        if self.is_bigwig:
            fetcher = bigWigIntervalQuery
        elif self.is_bigbed:
            fetcher = bigBedCoverageIntervals

        # check the coordinate inputs
        if bins < 0 and len(np.unique(ends_ - starts_)) != 1:
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
            chromSize = bbiChromSize(self.bbi, chromName)
            if chromSize == 0:
                raise KeyError("Chromosome not found: {}".format(chroms_[i]))
            start = starts_[i]
            end = ends_[i]
            if start > chromSize:
                raise ValueError(
                    "Start exceeds the chromosome length, {}.".format(chromSize))
            if bins < 1:
                array_query_full(
                    out[i, :], bins, self.bbi, fetcher,
                    chromName, start, end, chromSize, oob
                )
            else:
                try:
                    summary_type = BBI_SUMMARY_TYPES[summary]
                except KeyError:
                    raise ValueError(
                        'Invalid summary type "{}". Must be one of: {}.'.format(
                            summary,
                            set(BBI_SUMMARY_TYPES.keys()))
                        )
                array_query_summarized(
                    out[i, :], bins, self.bbi, fetcher,
                    chromName, start, end, chromSize, oob, summary_type,
                    exact
                )
        return out

    def fetch_intervals(self, str chrom, int start, int end, bint iterator=False):
        """
        Return an iterator or data frame of feature intervals overlapping a
        specified genomic query interval.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            Start coordinate.
        end : int
            End coordinate. If end is less than zero, the end is set to the
            chromosome size.
        iterator : bool, optional
            If True, return an iterator that provides records as tuples. For
            bigBeds, extra BED fields will be returned as unparsed strings. If
            False, return a DataFrame with all columns parsed according to the
            file's autosql schema.

        Returns
        -------
        pandas.DataFrame or BigWigIntervalIterator or BigBedIntervalIterator
        """
        if self.bbi == NULL:
            raise OSError("File closed")

        # find the chromosome
        cdef bytes chromName = chrom.encode('ascii')
        cdef int chromSize = bbiChromSize(self.bbi, chromName)
        if chromSize == 0:
            raise KeyError("Chromosome not found: {}".format(chrom))

        # check the coordinates
        if end < 0:
            end = chromSize
        if start > chromSize:
            raise ValueError(
                "Start exceeds the chromosome length, {}.".format(chromSize))
        cdef int length = end - start
        if length < 0:
            raise ValueError(
                "Interval cannot have negative length:"
                " start = {}, end = {}.".format(start, end))

        # clip the query range
        cdef int validStart = start
        cdef int validEnd = end
        if start < 0:
            validStart = 0
        if end > chromSize:
            validEnd = chromSize

        if self.is_bigwig:
            it = BigWigIntervalIterator(self, chromName, validStart, validEnd)
        else:
            it = BigBedIntervalIterator(self, chromName, validStart, validEnd)

        if iterator:
            return it
        else:
            try:
                import pandas as pd
            except ImportError:
                raise ImportError("fetch_intervals requires pandas")

            df = pd.DataFrame(list(it), columns=self.schema['columns'])
            for col, dtype in self.schema['dtypes'].items():
                try:
                    df.loc[:, col] = df[col].astype(dtype)
                except:
                    pass

        return df

    def fetch_summaries(self, str chrom, int start, int end, int zoom=0):
        """
        Return a DataFrame of complete summary statistic bins overlapping a
        genomic query interval at a specific zoom level.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            Start coordinate.
        end : int
            End coordinate. If end is less than zero, the end is set to the
            chromosome size.
        zoom : int, optional
            Zoom level. Default is 0.

        Returns
        -------
        pandas.DataFrame

        Notes
        -----
        The summary fields are 'validCount', 'min', 'max', 'sum', 'sumSquares'.

        See Also
        --------
        fetch : Fetch the signal track of a single interval
        fetch_intervals : Fetch feature intervals overlapping a genomic interval
        """
        if self.bbi == NULL:
            raise OSError("File closed")

        cdef bytes chromName = chrom.encode('ascii')

        if zoom < 0 or zoom >= len(self.zooms):
            raise ValueError(f"Zoom level {zoom} does not exist")

        cdef bbiZoomLevel *zoomObj = NULL
        cdef bbiZoomLevel *level = self.bbi.levelList
        cdef int i = 0
        while level != NULL:
            if zoom == i:
                zoomObj = level
                break
            level = level.next
            i += 1

        return get_zoom_table(self.bbi, zoomObj, chromName, start, end)

    def best_zoom(self, str chrom, int start, int end, int bins):
        cdef int baseSize = end - start
        cdef int stepSize = baseSize // bins
        cdef int zoomLevel = stepSize // 2
        if zoomLevel < 0:
            zoomLevel = 0

        cdef bbiZoomLevel *zoomObj = bbiBestZoom(self.bbi.levelList, zoomLevel)
        if zoomObj == NULL:
            return None

        cdef int binsize = zoomObj.reductionLevel
        cdef int i = 0
        for i, resolution in enumerate(self.zooms):
            if resolution == binsize:
                break
        return i


cdef class BigWigIntervalIterator:

    cdef str chrom
    cdef int valid_start
    cdef int valid_end
    cdef bbiInterval *interval
    cdef lm *lm

    def __init__(self, BBIFile fp, bytes chromName, int validStart, int validEnd):
        if fp.closed:
            raise OSError("File closed")
        self.chrom = chromName.decode('ascii')
        self.valid_start = validStart
        self.valid_end = validEnd

        # interval list is allocated out of lm
        self.lm = lmInit(0)
        self.interval = bigWigIntervalQuery(
            fp.bbi, chromName, validStart, validEnd, self.lm
        )

    def __iter__(self):
        return self

    def __next__(self):
        if self.interval == NULL:
            raise StopIteration

        cdef tuple out = (
            self.chrom, self.interval.start, self.interval.end, self.interval.val
        )
        self.interval = self.interval.next

        return out

    def __dealloc__(self):
        if self.interval != NULL:
            self.interval = NULL
        if self.lm != NULL:
            lmCleanup(&self.lm)


cdef class BigBedIntervalIterator:

    cdef str chrom
    cdef int valid_start
    cdef int valid_end
    cdef bigBedInterval *interval
    cdef lm *lm

    def __init__(self, BBIFile fp, bytes chromName, int validStart, int validEnd):
        if fp.closed:
            raise OSError("File closed")
        self.chrom = chromName.decode('ascii')
        self.valid_start = validStart
        self.valid_end = validEnd

        # interval list is allocated out of lm
        self.lm = lmInit(0)
        self.interval = bigBedIntervalQuery(
            fp.bbi, chromName, validStart, validEnd, 0, self.lm
        )

    def __iter__(self):
        return self

    def __next__(self):
        if self.interval == NULL:
            raise StopIteration

        cdef tuple restTup
        cdef str rest
        cdef char *cRest = self.interval.rest
        if cRest != NULL:
            rest = (<bytes>cRest).decode('ascii')
            restTup = tuple(rest.split('\t'))
        else:
            restTup = ()

        out = (self.chrom, self.interval.start, self.interval.end, *restTup)
        self.interval = self.interval.next
        return out

    def __dealloc__(self):
        if self.interval != NULL:
            self.interval = NULL
        if self.lm != NULL:
            lmCleanup(&self.lm)


cdef get_zoom_table(
    bbiFile *bbi,
    bbiZoomLevel *zoomObj,
    bytes chromName,
    int start,
    int end,
):
    """
    Get the summary elements overlapping the query interval for the given zoom level
    """
    cdef int chromId = bbiChromId(bbi, chromName)
    cdef int chromSize = bbiChromSize(bbi, chromName)
    if chromSize == 0:
        raise KeyError(f"Chromosome not found: {chromName}")

    # Clip the query range
    cdef int validStart = start
    cdef int validEnd = end
    if start < 0:
        validStart = 0
    if end > chromSize:
        validEnd = chromSize

    # Find appropriate zoom-level summary data and populate a list of summaries
    # summList is allocated and freed
    cdef bbiSummary *summList = bbiSummariesInRegion(
        zoomObj, bbi, chromId, validStart, validEnd
    )
    cdef int summ_start
    cdef int summ_end
    cdef int summ_validCount
    cdef float summ_min
    cdef float summ_max
    cdef float summ_sum
    cdef float summ_sumSquares

    cdef list summary_bins = []
    cdef bbiSummary *summ = summList
    while summ != NULL:
        summ_start = summ.start
        summ_end = summ.end
        summ_validCount = summ.validCount
        summ_min = summ.minVal
        summ_max = summ.maxVal
        summ_sum = summ.sumData
        summ_sumSquares = summ.sumSquares
        summary_bins.append(
            (
                chromName.decode('ascii'),
                summ_start,
                summ_end,
                summ_validCount,
                summ_min,
                summ_max,
                summ_sum,
                summ_sumSquares,
            )
        )
        summ = summ.next
    slFreeList(&summList)

    import pandas as pd

    return pd.DataFrame(
        summary_bins,
        columns=[
            "chrom",
            "start",
            "end",
            "validCount",
            "min",
            "max",
            "sum",
            "sumSquares",
        ],
    )


cdef inline void array_query_full(
    np.ndarray[np.double_t, ndim=1] out,
    int nbins,
    bbiFile *bbi,
    BbiFetchIntervals fetchIntervals,
    bytes chromName,
    int start,
    int end,
    int chromSize,
    double oob
):
    # Clip the query range
    cdef int validStart = start
    cdef int validEnd = end
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
    bbiSummaryType summaryType,
    boolean exact
):

    # Clip the query range
    cdef int validStart = start
    cdef int validEnd = end
    if start < 0:
        validStart = 0
    if end > chromSize:
        validEnd = chromSize

    # If not exact, get the closest zoom level less than what we're looking for
    cdef int baseSize
    cdef int stepSize
    cdef int zoomLevel
    cdef bbiZoomLevel *zoomObj = NULL
    baseSize = end - start
    stepSize = baseSize // nbins
    if not exact:
        zoomLevel = stepSize // 2
        if zoomLevel < 0:
            zoomLevel = 0
        zoomObj = bbiBestZoom(bbi.levelList, zoomLevel)

    # Allocate the elements array and create and populate summary elements
    cdef boolean result = False
    cdef bbiSummaryElement *elements = NULL
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
   int nbins
):
    # Look up region in index and get data at given zoom level.
    # Summarize this data in the summary array.

    cdef int chromId = bbiChromId(bbi, chromName)
    if chromId < 0:
        return False

    # Find appropriate zoom-level summary data
    # summList is allocated
    cdef bbiSummary *summList = bbiSummariesInRegion(
        zoom, bbi, chromId, validStart, validEnd
    )

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
    int nbins
):
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
