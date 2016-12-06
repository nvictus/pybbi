#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os.path as op

thisdir = op.dirname(op.realpath(__file__))


ext_modules = [
    Extension(
        name='kent.bbi', 
        sources=['kent/bbi.pyx'],
        library_dirs=[
            op.join(thisdir, 'src/x86_64'),
        ],
        libraries=[
            'kent',
        ],
        include_dirs=[
            np.get_include(),
            op.join(thisdir, 'include'),
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
    ext_modules=cythonize(ext_modules),
)
