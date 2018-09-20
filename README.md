# pybbi #

Python interface to Jim Kent's big binary file (bbi) \[[1](#ref1)\] library from the [UCSC Genome Browser source tree](https://github.com/ucscGenomeBrowser/kent) using Cython.

This provides read-level access to local and remote bigWig and bigBed files but no write capabilitites. The main feature is fast retrieval of range queries into numpy arrays.

## API ##

### Introspection

These accept a local file path or URL.

- `bbi.is_bbi(path)` --> `bool`
- `bbi.is_bigwig(path)` --> `bool`
- `bbi.is_bigbed(path)` --> `bool`
- `bbi.chromsizes(path)` --> `OrderedDict`
- `bbi.zooms(path)` --> `list`
- `bbi.info(path)` --> `dict`

### Signal query

These accept either a bigWig or bigBed file path / URL. The signal of a bigBed file is the coverage of its intervals.

- `bbi.fetch(path, chrom, start, end, [bins [, missing [, oob]]])` --> 1D numpy array
- `bbi.stackup(path, chroms, starts, ends, [bins [, missing [, oob]]])` --> 2D numpy array ("stacked heatmap")

Data "**summary**" querying is supported by specifying the number of `bins` for coarse graining. Currently, only the mean statistic is supported. **Missing** data can be filled with a custom fill value, `missing` (default = 0). Finally, **out-of-bounds** ranges (i.e. `start` less than zero or `end` greater than the chromosome length) are permitted because of their utility e.g., for generating vertical heatmap stacks centered at specific genomic features. A separate custom fill value, `oob` can be provided for out-of-bounds positions (default = NaN).

### Interval query

Accepts either a bigWig or bigBed file path / URL.

- `bbi.fetch_intervals(path, chrom, start, end)` --> iterator


See the docstrings for complete documentation.


## Installation ##

Requires
- Linux/MacOS
- C compiler, zlib, pthreads, libpng, openssl, make
- Python 2.7/3.3+
- `numpy` and `cython`

`pybbi` ships with (slightly modified) kent utils source, which it will compile before building the extension module.

```
$ pip install git+https://github.com/nvictus/pybbi.git
```

Or clone and install in development mode:

```
$ git clone https://github.com/nvictus/pybbi.git
$ cd pybbi
$ pip install -e .
```

## Related projects ##

- [libBigWig](https://github.com/dpryan79/libBigWig): Alternative C library for bigWig and bigBed files by Devon Ryan
- [pyBigWig](https://github.com/dpryan79/pyBigWig): Python bindings for `libBigWig` by the same author
- [bw-python](https://github.com/brentp/bw-python): Alternative Python wrapper to `libBigWig` by Brent Pederson
- [bx-python](https://github.com/bxlab/bx-python): Python bioinformatics library from James Taylor's group that includes tools for bbi files.

## References ##

<a id="ref1">[1]</a>: http://bioinformatics.oxfordjournals.org/content/26/17/2204.full

## Troubleshooting ##

On OSX, you may get errors about missing header files (e.g., `png.h`, `openssl/sha.h`), which even if installed may not be located in standard include locations. Either [create the required symlinks](https://www.anintegratedworld.com/mac-osx-fatal-error-opensslsha-h-file-not-found/) or update the `C_INCLUDE_PATH` environment variable accordingly before installing pybbi.

```bash
export C_INCLUDE_PATH="/usr/local/include/libpng:/usr/local/opt/openssl/include:$C_INCLUDE_PATH"
```

### Notes

[pyBigWig](https://github.com/dpryan79/pyBigWig) now also provides numpy-based retrieval and bigBed support.

Unfortunately, Kent's C source is not well-behaved library code, as it is littered with error calls that call `exit()`. I've added measures to `pybbi` to pre-empt common input errors, but if an internal error does get thrown, it will crash your interpreter instance. Check out [@dpryan79](https://github.com/dpryan79)'s [libBigWig](https://github.com/dpryan79/libBigWig) for an alternative and dedicated C library for big binary files.
