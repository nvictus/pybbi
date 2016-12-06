#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
from subprocess import check_call
import os.path as op
import numpy as np
import re
import io

thisdir = op.dirname(op.realpath(__file__))


class build_ext(_build_ext):
    def run(self):
        # First, compile our C library
        print("Compiling libkent...")
        check_call(['make', 'build-c'])
        # Now, proceed to build extension modules
        _build_ext.run(self)


ext_modules = [
    Extension(
        name='bbi.cbbi', 
        sources=[
            op.join(thisdir, 'bbi/cbbi.pyx')
        ],
        library_dirs=[
            op.join(thisdir, 'src/x86_64'),
        ],
        libraries=[
            'c', 'z', 'pthread', 'kent',
        ],
        include_dirs=[
            np.get_include(),
            op.join(thisdir, 'include'),
        ],
    ),
]


def read(*parts, **kwargs):
    encoding = kwargs.pop('encoding', 'utf-8')
    filepath = op.join(op.dirname(__file__), *parts)
    return io.open(filepath, encoding=encoding).read()


def get_version(pkg):
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        read(pkg, '__init__.py'),
        re.MULTILINE).group(1)
    return version


setup(
    name='pybbi',
    url='https://github.com/nvictus/pybbi',
    version=get_version('bbi'),
    packages=['bbi'],
    zip_safe=False,
    install_requires=['numpy', 'cython'],
    tests_require=['pytest'],
    ext_modules=cythonize(ext_modules),
    cmdclass={
        'build_ext': build_ext
    }
)
