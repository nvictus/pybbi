import numpy as np

from . import cbbi

__all__ = ["info", "chromsizes", "zooms", "fetch_intervals", "fetch", "stackup"]


def documented_by(original):
    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


def chromsizes(inFile):
    """
    Fetch the chromosome list of a bbi file. Returns an ordered dictionary of
    chromosome names mapped to their sizes in bp.

    Parameters
    ----------
    inFile : str
        Path to BigWig or BigBed file.

    Returns
    -------
    dict[str, int]

    """
    with cbbi.open(inFile) as f:
        return f.chromsizes


def zooms(inFile):
    """
    Fetch the zoom levels of a bbi file. Returns a list of "reduction levels",
    i.e. the number of bases per summary item, i.e. the bin size.

    Parameters
    ----------
    inFile : str
        Path to BigWig or BigBed file.

    Returns
    -------
    list[int]

    """
    with cbbi.open(inFile) as f:
        return f.zooms


def info(inFile):
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
    with cbbi.open(inFile) as f:
        return f.info


@documented_by(cbbi.BBIFile.fetch)
def fetch(
    inFile,
    chrom,
    start,
    end,
    bins=-1,
    missing=0.0,
    oob=np.nan,
    summary="mean",
    exact=False,
):
    with cbbi.open(inFile) as f:
        return f.fetch(chrom, start, end, bins, missing, oob, summary, exact)


@documented_by(cbbi.BBIFile.stackup)
def stackup(
    inFile,
    chroms,
    starts,
    ends,
    bins=-1,
    missing=0.0,
    oob=np.nan,
    summary="mean",
    exact=False,
):
    with cbbi.open(inFile) as f:
        return f.stackup(chroms, starts, ends, bins, missing, oob, summary, exact)


@documented_by(cbbi.BBIFile.fetch_intervals)
def fetch_intervals(inFile, chrom, start, end, iterator=True):
    with cbbi.open(inFile) as f:
        return f.fetch_intervals(chrom, start, end, iterator)


@documented_by(cbbi.BBIFile.fetch_summaries)
def fetch_summaries(inFile, chrom, start, end, zoom=0):
    with cbbi.open(inFile) as f:
        return f.fetch_summaries(chrom, start, end, zoom)
