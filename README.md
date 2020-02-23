# pybbi #

Python interface to Jim Kent's big binary file (bbi) \[[1](#ref1)\] library from the [UCSC Genome Browser source tree](https://github.com/ucscGenomeBrowser/kent) using Cython.

This provides read-level access to local and remote bigWig and bigBed files but no write capabilitites. The main feature is fast retrieval of range queries into numpy arrays.


## Installation ##

Wheels for `pybbi` are available on PyPI for Pythons 3.5, 3.6, 3.7, 3.8 on Linux and Mac OSX.

```
$ pip install pybbi
```

## API ##

### Introspection

These accept a local file path or URL.

- `bbi.is_bbi(path)` --> `bool`
- `bbi.is_bigwig(path)` --> `bool`
- `bbi.is_bigbed(path)` --> `bool`
- `bbi.chromsizes(path)` --> `OrderedDict`
- `bbi.zooms(path)` --> `list`
- `bbi.info(path)` --> `dict`

### Array output

These accept either a bigWig or bigBed file path / URL. The signal of a bigBed file is the genomic coverage of its intervals.

For a single range query:
- `bbi.fetch(path, chrom, start, end, [bins [, missing [, oob, [, summary]]]])` --> 1D numpy array

For a list of equal-length segments (i.e. to produce a stacked heatmap):
- `bbi.stackup(path, chroms, starts, ends, [bins [, missing [, oob, [, summary]]]])` --> 2D numpy array

**Summary** querying is supported by specifying the number of `bins` for coarsegraining. The `summary` statistic can be one of: 'mean', 'min', 'max', 'cov', or 'std'. (default = 'mean').

**Missing** data can be filled with a custom fill value, `missing` (default = 0). 

**Out-of-bounds** ranges (i.e. `start` less than zero or `end` greater than the chromosome length) are permitted because of their utility e.g., for generating vertical heatmap stacks centered at specific genomic features. A separate custom fill value, `oob` can be provided for out-of-bounds positions (default = NaN).

### Interval output

Accepts either a bigWig or bigBed file path / URL.

- `bbi.fetch_intervals(path, chrom, start, end)` --> iterator

See the docstrings for complete documentation.

## Related projects ##

- [libBigWig](https://github.com/dpryan79/libBigWig): Alternative C library for bigWig and bigBed files by Devon Ryan
- [pyBigWig](https://github.com/dpryan79/pyBigWig): Python bindings for `libBigWig` by the same author
- [bw-python](https://github.com/brentp/bw-python): Alternative Python wrapper to `libBigWig` by Brent Pederson
- [bx-python](https://github.com/bxlab/bx-python): Python bioinformatics library from James Taylor's group that includes tools for bbi files.

This library provides bindings to the reference UCSC bbi library code. Check out [@dpryan79](https://github.com/dpryan79)'s [libBigWig](https://github.com/dpryan79/libBigWig) for an alternative and dedicated C library for big binary files. pyBigWig also provides numpy-based retrieval and bigBed support.

## References ##

<a id="ref1">[1]</a>: http://bioinformatics.oxfordjournals.org/content/26/17/2204.full

## From source ##

If wheels for your platform or Python version aren't available or you want to develop, you'll need to install `pybbi` from source. The source distribution on PyPI ships with (slightly modified) kent utils source, which will compile before the extension module is built.

Requires
- Linux/MacOS
- C compiler, zlib, pthreads, libpng, openssl, make
- Python 2.7/3.4+
- `numpy` and `cython`

On fresh Ubuntu instance, you'll need `build-essential`, `make`, `zlib1g-dev`, `libssl-dev`, `libpng16-dev`. It seems to work on the Windows Subsystem for Linux too.

On a Centos/RedHat (rpm) system you'll need `gcc`, `make`, `zlib-devel`, `openssl-devel`, `libpng-devel`.

For development, clone the repo and install in editable mode:

```
$ git clone https://github.com/nvictus/pybbi.git
$ cd pybbi
$ pip install -e .
```

### Troubleshooting

On OSX, you may get errors about missing header files (e.g., `png.h`, `openssl/sha.h`), which even if installed may not be located in standard include locations. Either [create the required symlinks](https://www.anintegratedworld.com/mac-osx-fatal-error-opensslsha-h-file-not-found/) or update the `C_INCLUDE_PATH` environment variable accordingly before installing pybbi.

```bash
export C_INCLUDE_PATH="/usr/local/include/libpng:/usr/local/opt/openssl/include:$C_INCLUDE_PATH"
```

### Notes

Unfortunately, Kent's C source is not well-behaved library code, as it is littered with error calls that call `exit()`. `pybbi` will catch and pre-empt common input errors, but if somehow an internal error does get raised, it will terminate your interpreter instance.
