#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from subprocess import check_call
from distutils import log
import os.path as op
import glob
import sys
import os
import re
import io


thisdir = op.dirname(op.realpath(__file__))

# https://solitum.net/openssl-os-x-el-capitan-and-brew/
INCLUDE_DIRS = [
    op.join(sys.prefix, 'include'),
    op.join(sys.prefix, 'include', 'openssl'),
    '/usr/include',
    '/usr/include/openssl',
    '/usr/local/opt/openssl/lib',  # darwin
]
INCLUDE_DIRS += glob.glob(op.join(sys.prefix, 'include', 'libpng*'))
INCLUDE_DIRS += glob.glob('/usr/local/include/libpng*')  # darwin

LIBRARY_DIRS = [
    op.join(sys.prefix, 'lib'),
    '/lib',
    '/usr/lib',
    '/usr/local/lib'
    '/usr/local/opt/openssl/lib',  # darwin
]


include_dirs = []
for inc_dir in INCLUDE_DIRS:
    if op.isdir(inc_dir):
        include_dirs.append(inc_dir)
        break
pth = os.environ.get("CPATH")
if pth:
    include_dirs += pth.split(':')
os.environ["CPATH"] = ':'.join(include_dirs)


library_dirs = []
for lib_dir in LIBRARY_DIRS:
    if op.isdir(lib_dir):
        library_dirs.append(lib_dir)
        break
pth = os.environ.get("LIBRARY_PATH")
if pth:
    library_dirs += pth.split(':')
os.environ["LIBRARY_PATH"] = ':'.join(library_dirs)


os.environ["LDFLAGS"] = '-Wl,--no-as-needed ' + os.environ.get("LDFLAGS", "")


class lazylist(list):
    '''
    Used to delay the build-time Cython and numpy imports required to configure
    our list of extension modules until after setup_requires has been processed.

    '''
    def __init__(self, callback):
        self._list = None
        self.callback = callback

    def _cached_list(self):
        if self._list is None:
            self._list = self.callback()
        return self._list

    def __len__(self):
        return len(self._cached_list())

    def __iter__(self):
        for e in self._cached_list():
            yield e

    def __getitem__(self, i):
        return self._cached_list()[i]


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


class build_ext(_build_ext):
    def run(self):
        os.environ['SETUP_PY'] = '1'

        # First, compile our C library: libkent.a
        import sysconfig
        log.info(sysconfig.get_config_vars())
        log.info("CPATH: " + os.environ.get("CPATH", ""))
        log.info("LIBRARY_PATH: " + os.environ.get("LIBRARY_PATH", ""))
        log.info("Compiling libkent archive...")
        check_call(['make', 'build-c'])

        # Now, proceed to build extension modules
        log.info("Building extension module...")
        _build_ext.run(self)


def get_ext_modules():
    from Cython.Build import cythonize
    import numpy

    ext_modules = [
        Extension(
            name='bbi.cbbi',
            sources=[
                op.join(thisdir, 'bbi/cbbi.pyx')
            ],
            libraries=[
                'c', 'z', 'pthread', 'ssl', 'crypto', 'png', 'kent'
            ],
            library_dirs=[
                op.join(thisdir, 'src/x86_64'),
            ] + library_dirs,
            include_dirs=[
                numpy.get_include(),
                op.join(thisdir, 'include'),
            ] + include_dirs,
        ),
    ]
    return cythonize(ext_modules)


setup(
    name='pybbi',
    author='Nezar Abdennur',
    author_email='nabdennur@gmail.com',
    license='MIT',
    version=get_version('bbi'),
    packages=['bbi'],
    description='Python bindings to UCSC Big Binary (bigWig/bigBed) file library',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url='https://github.com/nvictus/pybbi',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    zip_safe=False,
    setup_requires=[
        'setuptools>=18.0',
        'cython',
        'numpy',
    ],
    install_requires=[
        'six',
        'numpy'
    ],
    tests_require=[
        'pytest'
    ],
    ext_modules=lazylist(get_ext_modules),
    cmdclass={
        'build_ext': build_ext
    }
)
