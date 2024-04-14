# pybbi #

![Build Status](https://github.com/nvictus/pybbi/actions/workflows/ci.yml/badge.svg)
[![DOI](https://zenodo.org/badge/58960207.svg)](https://zenodo.org/doi/10.5281/zenodo.10382980)

Python interface to Jim Kent's Big Binary Indexed file (BBI) \[[1](#ref1)\] library from the [UCSC Genome Browser source tree](https://github.com/ucscGenomeBrowser/kent) using Cython.

This provides read-level access to local and remote bigWig and bigBed files but no write capabilitites. The main feature is fast retrieval of range queries into numpy arrays.


## Installation ##

Wheels for `pybbi` are available on PyPI for Python 3.8, 3.9, 3.10, 3.11 on Linux (x86_64 and aarch64) and Mac OSX (x86_64/Intel). Apple Silicon (arm64) wheels will be made available once M1 runners are available in GitHub Actions.

```
$ pip install pybbi
```

## API ##

The `bbi.open` function returns a `BBIFile` object.

```
bbi.open(path) -> BBIFile
```

`path` can be a local file path (bigWig or bigBed) or a URL. `BBIFile` objects are context managers and can be used in a `with` statement to clean up resources without calling `BBIFile.close()`.

```python
>>> with bbi.open('bigWigExample.bw') as f:
...     x = f.fetch('chr21', 1000000, 2000000, bins=40)
```

### Introspection

```
BBIFile.is_bigwig -> bool
BBIFile.is_bigbed -> bool
BBIFile.chromsizes -> OrderedDict
BBIFile.zooms -> list
BBIFile.info -> dict
BBIFile.schema -> dict
BBIFile.read_autosql() -> str
```

Note: `BBIFile.schema['dtypes']` provides numpy data types for the fields in a bigWig or bigBed (matched from the autoSql definition).


### Interval output

The actual interval records in a bigWig or bigBed can be retrieved as a pandas dataframe or as an iterator over records as tuples. The pandas output is parsed according to the file's schema.

```
BBIFile.fetch_intervals(chrom, start, end) -> pandas.DataFrame
BBIFile.fetch_intervals(chrom, start, end, iterator=True) -> interval iterator
```

Summary bin records at each zoom level are also accessible.

```
BBIFile.fetch_summaries(chrom, start, end, zoom) -> pandas.DataFrame
```

### Array output

Retrieve quantitative signal as an array. The signal of a bigWig file is obtained from its "value" field. The signal of a bigBed file is obtained from the genomic coverage of its intervals.

For a single range query:
```
BBIFile.fetch(chrom, start, end, [bins [, missing [, oob, [, summary]]]]) -> 1D numpy array
```

To produce a stacked heatmap from a list of (1) equal-length intervals or (2) arbitrary-length intervals with `bins` specified:
```
BBIFile.stackup(chroms, starts, ends, [bins [, missing [, oob, [, summary]]]]) -> 2D numpy array
```

* **Summary** querying is supported by specifying the number of `bins` for coarsening. The `summary` statistic can be one of: 'mean', 'min', 'max', 'cov', 'std', 'or 'sum'. (default = 'mean'). Intervals need not have the same length, in which case the data from each interval will be interpolated to the same number of bins (e.g., gene bodies).

* **Missing** data can be filled with a custom fill value, `missing` (default = 0). 

* **Out-of-bounds** ranges (i.e. `start` less than zero or `end` greater than the chromosome length) are permitted because of their utility e.g., for generating vertical heatmap stacks centered at specific genomic features. A separate custom fill value, `oob` can be provided for out-of-bounds positions (default = NaN).

### Function API

The original function-based API is still available:

```python
bbi.is_bbi(path: str) -> bool
bbi.is_bigwig(path: str) -> bool
bbi.is_bigbed(path:str) -> bool
bbi.chromsizes(path: str) -> OrderedDict
bbi.zooms(path: str) -> list
bbi.info(path: str) -> dict
bbi.fetch_intervals(path: str, chrom: str, start: int, end: int, iterator: bool) -> Union[Iterable, pd.DataFrame]
bbi.fetch(path: str, chrom: str, start: int, end: int, [bins: int [, missing: float [, oob: float, [, summary: str]]]]) -> np.array[1, 'float64']
bbi.stackup(path: str, chroms: np.array, starts: np.array, ends: np.array, [bins: int [, missing: float [, oob: float, [, summary: str]]]]) -> np.array[2, 'float64']
```

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
- Platform: Linux or Darwin (Windows Subsystem for Linux seems to work too)
- pthreads, zlib, libpng, openssl, make, pkg-config
- Python 3.6+
- `numpy` and `cython`

For example, on a fresh Ubuntu instance, you'll need `build-essential`, `make`, `pkg-config`, `zlib1g-dev`, `libssl-dev`, `libpng16-dev`.

On a Centos/RedHat (rpm) system you'll need `gcc`, `make`, `pkg-config`, `zlib-devel`, `openssl-devel`, `libpng-devel`.

On a Mac, you'll need Xcode and to `brew install pkg-config openssl libpng`.

For development, clone the repo and install in editable mode:

```
$ git clone https://github.com/nvictus/pybbi.git
$ cd pybbi
$ pip install -e .
```

You can use the `ARCH` environment variable to specify a target architecture or `ARCHFLAGS` on a Mac.

### Notes

Unfortunately, Kent's C source is not well-behaved library code, as it is littered with error calls that call `exit()`. `pybbi` will catch and pre-empt common input errors, but if somehow an internal error does get raised, it will terminate your interpreter instance.
