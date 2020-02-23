#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from subprocess import check_call
from distutils import log
import os.path as op
import sys
import os
import re
import io


thisdir = op.dirname(op.realpath(__file__))


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
        # import sysconfig
        # log.info(sysconfig.get_config_vars())
        # log.info("CPATH: " + os.environ.get("CPATH", ""))
        # log.info("LIBRARY_PATH: " + os.environ.get("LIBRARY_PATH", ""))
        log.info("Compiling libkent archive...")
        check_call(['make', 'build-ucsc'])

        # Now, proceed to build extension modules
        log.info("Building extension module...")
        _build_ext.run(self)


def get_ext_modules():
    from Cython.Build import cythonize
    import numpy
    import pkgconfig
    import sysconfig

    # https://solitum.net/openssl-os-x-el-capitan-and-brew/
    if sys.platform == "darwin":
        s = '/usr/local/opt/openssl/lib/pkgconfig'
        old = os.environ.get('PKG_CONFIG_PATH')
        if old:
            s = old + ':' + s
        os.environ['PKG_CONFIG_PATH'] = s
        os.environ['MACOSX_DEPLOYMENT_TARGET'] = \
            sysconfig.get_config_var('MACOSX_DEPLOYMENT_TARGET')
        os.environ['BLDSHARED'] = \
            'gcc -bundle -undefined dynamic_lookup -arch x86_64 -g'
        os.environ['LDSHARED'] = \
            'gcc -bundle -undefined dynamic_lookup -arch x86_64 -g'

    if sys.platform == "linux":
        s = '-Wl,--no-as-needed'
        old = os.environ.get("LDFLAGS")
        if old:
            s = s + ' ' + old
        os.environ["LDFLAGS"] = s

    d = pkgconfig.parse('zlib openssl libpng')

    ext_modules = [
        Extension(
            name='bbi.cbbi',
            sources=[
                op.join(thisdir, 'bbi/cbbi.pyx')
            ],
            libraries=[
                'kent',
            ] + d.pop('libraries', []),
            library_dirs=[
                op.join(thisdir, 'src/x86_64'),
            ] + d.pop('library_dirs', []),
            include_dirs=[
                numpy.get_include(),
                op.join(thisdir, 'include'),
            ] + d.pop('include_dirs', []),
            **d
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
        'pkgconfig'
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
