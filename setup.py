#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from subprocess import check_call
import os.path as op
import sys
import os
import re
import io


thisdir = op.dirname(op.realpath(__file__))


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


class lazy_list(list):
    '''
    Used to delay cython import until all setup_requires have been installed.

    '''
    def __init__(self, callback):
        self._list, self.callback = None, callback
    def c_list(self):
        if self._list is None:
            self._list = self.callback()
        return self._list
    def __iter__(self):
        for e in self.c_list():
            yield e
    def __getitem__(self, i):
        return self.c_list()[i]
    def __len__(self):
        return len(self.c_list())


# Paths to propagate to make to build libkent
extra_library_dirs = []
extra_include_dirs = []
for dirname in ['libpng16', 'openssl']:
    inc_dir = op.join(sys.prefix, 'include', dirname)
    if op.isdir(inc_dir):
        extra_include_dirs.append(inc_dir)


if sys.platform == "darwin":
    # https://solitum.net/openssl-os-x-el-capitan-and-brew/
    extra_library_dirs += [
        '/usr/local/opt/openssl/lib'
    ]
    extra_include_dirs += [
        '/usr/local/opt/openssl/include', 
        '/usr/local/include/libpng16'
    ]


class build_ext(_build_ext):
    def run(self):
        # First, compile our C library
        for lib_dir in extra_library_dirs[::-1]:
            os.environ["LIBRARY_PATH"] = lib_dir + ':' + os.environ.get("LIBRARY_PATH", "")
        for inc_dir in extra_include_dirs[::-1]:
            os.environ["C_INCLUDE_PATH"] = inc_dir + ':' + os.environ.get("C_INCLUDE_PATH", "")
        print("Compiling libkent...", file=sys.stderr)
        print("LIBRARY_PATH: " + os.environ.get("LIBRARY_PATH", ""), file=sys.stderr)
        print("C_INCLUDE_PATH: " + os.environ.get("C_INCLUDE_PATH", ""), file=sys.stderr)
        check_call(['make', 'build-c'])
        # Now, proceed to build extension modules
        _build_ext.run(self)

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Fix to work with bootstrapped numpy installation
        # http://stackoverflow.com/a/21621689/579416
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


def get_ext_modules():
    from Cython.Build import cythonize
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
                #np.get_include(),
                op.join(thisdir, 'include'),
            ] + extra_include_dirs,
        ),
    ]
    return cythonize(ext_modules)


setup(
    name='pybbi',
    url='https://github.com/nvictus/pybbi',
    version=get_version('bbi'),
    packages=['bbi'],
    zip_safe=False,
    setup_requires=[
        'setuptools>=18.0',
        'cython',
        'numpy',
    ],
    install_requires=['six', 'numpy'],
    tests_require=['pytest'],
    ext_modules=lazy_list(get_ext_modules),
    cmdclass={
        'build_ext': build_ext
    }
)
