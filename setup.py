import os
import pathlib
import subprocess
import sys

from setuptools import Extension, setup

# Determine architecture for build; defaults to machine architecture.
# cibuildwheel sets ARCHFLAGS on macos runners
# user can also set ARCH to override
MACHINE_ARCH = os.uname().machine
if sys.platform == "darwin" and "ARCHFLAGS" in os.environ:
    ARCH = os.environ["ARCHFLAGS"].split()[-1]
elif "ARCH" in os.environ:
    ARCH = os.environ["ARCH"]
else:
    ARCH = MACHINE_ARCH


if "sdist" in sys.argv:
    # Skip compilation when building a source distribution
    setup()

else:
    import numpy
    import pkgconfig
    from Cython.Build import cythonize
    from Cython.Distutils.build_ext import new_build_ext as cython_build_ext

    # 1. Compile the UCSC library (libkent.a)
    # Pass the target architecture to the makefile
    os.environ["MACHTYPE"] = ARCH

    # Platform and architecture-specific flags
    if sys.platform == "darwin":
        brew_optdir = None
        if ARCH.startswith("arm"):
            brew_optdir = "/opt/homebrew/opt"
        else:
            brew_optdir = "/usr/local/opt"
        os.environ["CFLAGS"] = f"-arch {ARCH}"
        if brew_optdir is not None:
            os.environ["LDFLAGS"] = f"-L{brew_optdir}/openssl/lib {os.environ.get('LDFLAGS', '')}"
            os.environ["C_INCLUDE_PATH"] = f"{brew_optdir}/openssl/include:{brew_optdir}/libpng/include:{os.environ.get('C_INCLUDE_PATH', '')}"
            os.environ["PKG_CONFIG_PATH"] = f"{brew_optdir}/openssl/lib/pkgconfig:{os.environ.get('PKG_CONFIG_PATH', '')}"

    # Parse pkg-config dependencies
    # This will let us know if something is missing
    deps = pkgconfig.parse("zlib openssl libpng")

    # Build the UCSC library
    if not pathlib.Path("src", ARCH, "libkent.a").exists():
        subprocess.run(["make", "clean-ucsc"])
        ret = subprocess.run(["make", "build-ucsc"], check=True)
        if ret.returncode != 0:
            raise RuntimeError("Failed to build UCSC library.")


    # 2. Compile cython extension module, link to libkent and other dependencies
    # Platform and architecture-specific linker flags
    if sys.platform == "darwin":
        os.environ["BLDSHARED"] = f'gcc -bundle -undefined dynamic_lookup -arch {ARCH} -g'
        os.environ["LDSHARED"] = f'gcc -bundle -undefined dynamic_lookup -arch {ARCH} -g'
    elif sys.platform == "linux":
        os.environ["LDFLAGS"] = f"-Wl,--no-as-needed {os.environ.get('LDFLAGS', '')}"
    
    # Configure and cythonize pyx to C extension module
    ext_modules = cythonize([
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
    ])


    setup(
        ext_modules=ext_modules,
        build_ext=cython_build_ext,
    )
