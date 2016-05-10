#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np
import os.path as op

thisdir = op.dirname(op.realpath(__file__))


ext_modules = [
    Extension(
        name='kent.bbi', 
        sources=['kent/bbi.pyx'],
        library_dirs=[
            op.join(thisdir, 'lib/x86_64'),
        ],
        libraries=[
            'kent',
        ],
        include_dirs=[
            np.get_include(),
            op.join(thisdir, 'include'),
            op.join(thisdir, 'src'),
        ],
    ),
]


setup(
    name='pykent',
    version='0.1',
    description='',
    author='N A',
    #author_email='',
    #url='',
    packages=['kent'],
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
)

