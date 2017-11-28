#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
from subprocess import check_call
import os.path as op
import numpy as np
import sys
import os
import re
import io

thisdir = op.dirname(op.realpath(__file__))


class build_ext(_build_ext):
    def run(self):
        # First, compile our C library
        for lib_dir in extra_library_dirs:
            os.environ["LIBRARY_PATH"] = lib_dir + ':' + os.environ.get("LIBRARY_PATH", "")
        for inc_dir in extra_include_dirs:
            os.environ["C_INCLUDE_PATH"] = inc_dir + ':' + os.environ.get("C_INCLUDE_PATH", "")
        print("Compiling libkent...", file=sys.stderr)
        print("LIBRARY_PATH: " + os.environ.get("LIBRARY_PATH", ""), file=sys.stderr)
        print("C_INCLUDE_PATH: " + os.environ.get("C_INCLUDE_PATH", ""), file=sys.stderr)
        check_call(['make', 'build-c'])

        # Now, proceed to build extension modules
        _build_ext.run(self)


extra_library_dirs = []
extra_include_dirs = []
if sys.platform == "darwin":
    # https://solitum.net/openssl-os-x-el-capitan-and-brew/
    extra_library_dirs += ['/usr/local/opt/openssl/lib']
    extra_include_dirs += ['/usr/local/opt/openssl/include']
elif sys.platform.startswith("linux"):
    for dirname in ['libpng16', 'openssl']:
        inc_dir = op.join(sys.prefix, 'include', dirname)
        if op.isdir(inc_dir):
            extra_include_dirs.append(inc_dir)


ext_modules = [
    Extension(
        name='bbi.cbbi', 
        sources=[
            op.join(thisdir, 'bbi/cbbi.pyx')
        ],
        library_dirs=[
            op.join(thisdir, 'src/x86_64'),
        ] + extra_library_dirs,
        libraries=[
            'c', 'z', 'pthread', 'kent',
        ],
        include_dirs=[
            np.get_include(),
            op.join(thisdir, 'include'),
        ] + extra_include_dirs,
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
