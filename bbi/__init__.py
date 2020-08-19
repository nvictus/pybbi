from .cbbi import (
    is_bbi,
    is_bigwig,
    is_bigbed,
    open,
    info,
    chromsizes,
    zooms,
    fetch,
    fetch_intervals,
    stackup
)

__all__ = [
    'is_bbi',
    'is_bigwig',
    'is_bigbed',
    'open',
    'info',
    'chromsizes',
    'zooms',
    'fetch',
    'fetch_intervals',
    'stackup',
]

__version__ = '0.2.3'
