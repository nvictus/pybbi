#-*- coding: utf-8 -*-
from .cbbi import (
    is_bbi,
    is_bigwig,
    is_bigbed,
    info,
    zooms,
    chromsizes,
    fetch_intervals,
    fetch,
    stackup,
    open
)

__all__ = [
    'is_bbi',
    'is_bigwig',
    'is_bigbed',
    'info',
    'zooms',
    'chromsizes',
    'fetch_intervals',
    'fetch',
    'stackup',
    'open'
]

__version__ = '0.2.1'
