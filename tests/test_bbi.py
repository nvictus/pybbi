import os.path as op
import numpy as np
import os

import bbi
import pytest

thisdir = op.dirname(op.realpath(__file__))
BW_FILE = op.join(thisdir, 'bigWigExample.bw')
BB_FILE = op.join(thisdir, 'bigBedExample.bb')
BW_URL = 'http://genome.ucsc.edu/goldenPath/help/examples/bigWigExample.bw'
BB_URL = 'http://genome.ucsc.edu/goldenPath/help/examples/bigBedExample.bb'

bbi_paths_and_urls = [BW_FILE, BW_URL, BB_FILE, BB_URL]
bbi_paths = [BW_FILE, BB_FILE]
bbi_urls = [BW_URL, BB_URL]


def test_sigs():
    assert bbi.is_bbi(BW_FILE)
    assert bbi.is_bigwig(BW_FILE)
    assert not bbi.is_bigbed(BW_FILE)
    with bbi.open(BW_FILE) as f:
        assert f.is_bigwig
        assert not f.is_bigbed
        assert not f.is_remote

    assert bbi.is_bbi(BW_URL)
    assert bbi.is_bigwig(BW_URL)
    assert not bbi.is_bigbed(BW_URL)
    with bbi.open(BW_URL) as f:
        assert f.is_bigwig
        assert not f.is_bigbed
        assert f.is_remote

    assert bbi.is_bbi(BB_FILE)
    assert not bbi.is_bigwig(BB_FILE)
    assert bbi.is_bigbed(BB_FILE)
    with bbi.open(BB_FILE) as f:
        assert not f.is_bigwig
        assert f.is_bigbed
        assert not f.is_remote

    assert bbi.is_bbi(BB_URL)
    assert not bbi.is_bigwig(BB_URL)
    assert bbi.is_bigbed(BB_URL)
    with bbi.open(BB_URL) as f:
        assert not f.is_bigwig
        assert f.is_bigbed
        assert f.is_remote


@pytest.mark.skipif(os.getenv("GITHUB_ACTIONS") == "true", reason="Skipped on GitHub Actions")
def test_aws_403_redirect():
    # See https://stat.ethz.ch/pipermail/bioc-devel/2016-May/009241.html
    url = 'https://www.encodeproject.org/files/ENCFF620UMO/@@download/ENCFF620UMO.bigWig'
    bbi.open(url).fetch('chr21', 0, 1000)


@pytest.mark.parametrize('uri', bbi_paths_and_urls)
def test_chromsizes(uri):
    with bbi.open(uri) as f:
        chromsizes = f.chromsizes
    assert len(chromsizes) == 1 and 'chr21' in chromsizes


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch(path):
    with bbi.open(path) as f:
        x = f.fetch('chr21', 0, 1000)
    assert len(x) == 1000

    with bbi.open(path) as f:
        x = f.fetch('chr21', 0, 1000, bins=10)
    assert len(x) == 10

    with pytest.raises(KeyError):
        with bbi.open(path) as f:
            f.fetch('chr1', 0, 1000)


def test_fetch_remote():
    x_local = bbi.open(BW_FILE).fetch('chr21', 0, 100)
    x_remote = bbi.open(BW_URL).fetch('chr21', 0, 100)
    assert np.allclose(x_local, x_remote, equal_nan=True)

    x_local = bbi.open(BB_FILE).fetch('chr21', 0, 100)
    x_remote = bbi.open(BB_URL).fetch('chr21', 0, 100)
    assert np.allclose(x_local, x_remote, equal_nan=True)


def test_fetch_remote_https():
    x_local = bbi.open(BW_FILE).fetch('chr21', 0, 100)
    x_remote = bbi.open(BW_URL.replace('http://', 'https://')).fetch(
        'chr21', 0, 100
    )
    assert np.allclose(x_local, x_remote, equal_nan=True)

    x_local = bbi.open(BB_FILE).fetch('chr21', 0, 100)
    x_remote = bbi.open(BB_URL.replace('http://', 'https://')).fetch(
        'chr21', 0, 100
    )
    assert np.allclose(x_local, x_remote, equal_nan=True)


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_missing(path):
    f = bbi.open(path)

    x = f.fetch('chr21', 0, 1000, oob=0)
    assert np.all(x[:10] == 0)
    x = f.fetch('chr21', 0, 1000, missing=np.nan)
    assert np.all(np.isnan(x[:10]))


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_oob(path):
    f = bbi.open(path)

    x = f.fetch('chr21', -10, 1000, oob=np.nan)
    assert np.all(np.isnan(x[:10]))
    x = f.fetch('chr21', -10, 1000, oob=0)
    assert np.all(x[:10] == 0)

    n = f.chromsizes['chr21']
    x = f.fetch('chr21', n - 1000, n + 10, oob=np.nan)
    assert np.all(np.isnan(x[-10:]))
    x = f.fetch('chr21', n - 1000, n + 10, oob=0)
    assert np.all(x[-10:] == 0)


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_intervals(path):
    it = bbi.fetch_intervals(path, 'chr21', 0, 1000)  # unmappable region
    assert len(list(it)) == 0

    it = bbi.fetch_intervals(path, 'chr21', 0, 10000000)
    assert len(list(it)) > 0

    with bbi.open(path) as f:
        df = f.fetch_intervals('chr21', 0, 1000)  # unmappable region
        assert len(df) == 0

    with bbi.open(path) as f:
        df = f.fetch_intervals('chr21', 0, 10000000)
        assert len(df) > 0


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_summary_stats(path):
    f = bbi.open(path)

    x = f.fetch('chr21', 20000000, 20001000, bins=10, summary='mean')
    y = f.fetch('chr21', 20000000, 20001000, bins=10)
    assert np.allclose(x, y)

    values = f.fetch('chr21', 20000000, 20001000)
    vmin = f.fetch(
        'chr21', 20000000, 20001000, bins=10, summary='min'
    ).min()
    assert np.isclose(vmin, np.min(values))
    vmax = f.fetch(
        'chr21', 20000000, 20001000, bins=10, summary='max'
    ).max()
    assert np.isclose(vmax, np.max(values))
    vsum = f.fetch('chr21', 20000000, 20001000, bins=100, summary='sum')
    values_sum_every_ten = np.reshape(values, (-1, 10)).sum(axis=-1)
    assert len(vsum) == len(values_sum_every_ten)
    assert np.allclose(vsum, values_sum_every_ten)

    with pytest.raises(ValueError):
        f.fetch('chr21', 20000000, 20001000, bins=10, summary='foo')


@pytest.mark.parametrize('path', bbi_paths)
def test_stackup(path):
    f = bbi.open(path)

    x = f.stackup(['chr21', 'chr21'], [0, 2000], [1000, 3000])
    assert x.shape == (2, 1000)

    x = f.stackup(['chr21', 'chr21'], [0, 2000], [1000, 3000], bins=10)
    assert x.shape == (2, 10)

    # unequal interval lengths
    with pytest.raises(ValueError):
        f.stackup(['chr21', 'chr21'], [0, 2000], [1000, 3500])

    x = f.stackup(['chr21', 'chr21'], [0, 2000], [1000, 3500], bins=10)
    assert x.shape == (2, 10)


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_zoom_records(path):
    with bbi.open(path) as f:
        df = f.fetch_summaries('chr21', 0, 10000000)
        assert len(df) > 0
