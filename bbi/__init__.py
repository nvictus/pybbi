#-*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
from .cbbi import (
    is_bbi,
    is_bigwig,
    is_bigbed,
    info,
    zooms,
    chromsizes,
    fetch,
    stackup,
    # fetch_intervals,
    BbiFile
)

__all__ = ['is_bbi', 'is_bigwig', 'is_bigbed', 'info', 'zooms', 'chromsizes',
           'fetch', 'stackup', 'BbiFile']  #  'bool' is not a type identifier

__version__ = '0.2.1-dev'
