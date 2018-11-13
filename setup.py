#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from subprocess import check_call, PIPE
import os.path as op
import subprocess
import sysconfig
import sys
import os
import re
import io


thisdir = op.dirname(op.realpath(__file__))


class delayedlist(list):
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


# Configuration for the static C library libkent.a
CONFIG_VARS = sysconfig.get_config_vars()
LIBDIR = CONFIG_VARS.get('LIBDIR', op.join(sys.prefix, 'lib'))
INCLUDEDIR = CONFIG_VARS.get('INCLUDEDIR', op.join(sys.prefix, 'include'))
MACHTYPE = subprocess.Popen(['uname', '-m'], stdout=PIPE).stdout.read()
if MACHTYPE:
    if sys.version_info[0] > 2:
        MACHTYPE = MACHTYPE.decode('utf-8')
    MACHTYPE = MACHTYPE.strip()

library_dirs = [
    LIBDIR,
    op.join(thisdir, 'src', MACHTYPE),
]

include_dirs = [
    INCLUDEDIR,
    op.join(thisdir, 'include'),
    op.join(thisdir, 'src'),
]

for dirname in ['libpng16', 'openssl']:
    inc_subdir = op.join(INCLUDEDIR, dirname)
    if op.isdir(INCLUDEDIR):
        include_dirs.append(inc_subdir)

if sys.platform == "darwin":
    # https://solitum.net/openssl-os-x-el-capitan-and-brew/
    library_dirs += [
        '/usr/local/opt/openssl/lib'
    ]
    include_dirs += [
        '/usr/local/opt/openssl/include',
        '/usr/local/include/libpng16'
    ]

_path = os.environ.get("LIBRARY_PATH", "")
if _path:
    _path = ':'.join(library_dirs + [_path])
else:
    _path = ':'.join(library_dirs)
os.environ["LIBRARY_PATH"] = _path

_path = os.environ.get("C_INCLUDE_PATH", "")
if _path:
    _path = ':'.join(include_dirs + [_path])
else:
    _path = ':'.join(include_dirs)
os.environ["C_INCLUDE_PATH"] = _path

os.environ["INC"] = ' '.join('-I' + d for d in include_dirs)

os.environ['LDFLAGS'] = '-lz -lc -lpthread'


# Configuration for the extension module
def get_ext_modules():
    from Cython.Build import cythonize
    import numpy

    ext_modules = [
        Extension(
            name='bbi.cbbi',
            sources=[op.join(thisdir, 'bbi', 'cbbi.pyx')],
            libraries=[
                'c',
                'z',
                'pthread',
                'kent',
            ],
            library_dirs=library_dirs,
            include_dirs=include_dirs + [numpy.get_include()],
        ),
    ]
    return cythonize(ext_modules)


class build_ext(_build_ext):
    def run(self):
        # First, build the static C library
        print("Compiling libkent...", file=sys.stderr)
        for var in ('LIBRARY_PATH', 'C_INCLUDE_PATH', 'LDFLAGS', 'INC'):
            print(var + ": "  + os.environ.get(var), file=sys.stderr)
        check_call(['make', 'build-c'])

        # Now, proceed to build the C extension module
        _build_ext.run(self)


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
    ext_modules=delayedlist(get_ext_modules),
    cmdclass={
        'build_ext': build_ext
    }
)
