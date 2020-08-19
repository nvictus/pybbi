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
    fetch,
    fetch_intervals,
    stackup,
)

del cbbi, _bbi

__version__ = '0.2.3'
