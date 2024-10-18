from .cbbi import (
    is_bbi,
    is_bigwig,
    is_bigbed,
    open,
)

from ._bbi import (
    info,
    chromsizes,
    zooms,
    fetch_intervals,
    fetch,
    stackup,
)

del cbbi, _bbi

__version__ = '0.4.1'
__all__ = [
    'is_bbi',
    'is_bigwig',
    'is_bigbed',
    'open',
    'info',
    'chromsizes',
    'zooms',
    'fetch_intervals',
    'fetch',
    'stackup',
]
