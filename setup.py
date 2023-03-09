import os
import os.path
import pathlib
import sys

import numpy
import pkgconfig
from Cython.Build import cythonize
from Cython.Distutils.build_ext import new_build_ext as cython_build_ext
from setuptools import Extension, setup

ARCH = os.uname().machine

deps = {}
if "sdist" not in sys.argv:
    deps = pkgconfig.parse("zlib openssl libpng")

    if not pathlib.Path("src", ARCH, "libkent.a").exists():
        raise RuntimeError(
            f"src/{ARCH}/libkent.a not found. "
            "Please run `make build-ucsc`."
        )


def get_extension_modules():
    ext_modules = [
        Extension(
            name="bbi.cbbi",
            sources=[os.path.join("bbi", "cbbi.pyx")],
            libraries=["kent"] + deps.pop("libraries", []),
            library_dirs=[
                os.path.join("src", ARCH),
            ]
            + deps.pop("library_dirs", []),
            include_dirs=[
                numpy.get_include(),
                "include",
            ]
            + deps.pop("include_dirs", []),
            **deps,
        ),
    ]
    if "sdist" in sys.argv:
        return ext_modules
    else:
        return cythonize(ext_modules)


setup(
    ext_modules=get_extension_modules(),
    build_ext=cython_build_ext,
)
