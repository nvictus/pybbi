# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import numpy as np
import requests

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

    assert bbi.is_bbi(BW_URL)
    assert bbi.is_bigwig(BW_URL)
    assert not bbi.is_bigbed(BW_URL)

    assert bbi.is_bbi(BB_FILE)
    assert not bbi.is_bigwig(BB_FILE)
    assert bbi.is_bigbed(BB_FILE)

    assert bbi.is_bbi(BB_URL)
    assert not bbi.is_bigwig(BB_URL)
    assert bbi.is_bigbed(BB_URL)


@pytest.mark.parametrize('uri', bbi_paths_and_urls)
def test_chromsizes(uri):
    chromsizes = bbi.chromsizes(uri)
    assert len(chromsizes) == 1 and 'chr21' in chromsizes


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch(path):
    x = bbi.fetch(path, 'chr21', 0, 1000)
    assert len(x) == 1000

    x = bbi.fetch(path, 'chr21', 0, 1000, bins=10)
    assert len(x) == 10

    with pytest.raises(KeyError):
        bbi.fetch(path, 'chr1', 0, 1000)


def test_fetch_remote():
    x_local = bbi.fetch(BW_FILE, 'chr21', 0, 100)
    x_remote = bbi.fetch(BW_URL, 'chr21', 0, 100)
    assert np.allclose(x_local, x_remote, equal_nan=True)

    x_local = bbi.fetch(BB_FILE, 'chr21', 0, 100)
    x_remote = bbi.fetch(BB_URL, 'chr21', 0, 100)
    assert np.allclose(x_local, x_remote, equal_nan=True)


def test_fetch_remote_https():
    x_local = bbi.fetch(BW_FILE, 'chr21', 0, 100)
    x_remote = bbi.fetch(
        BW_URL.replace('http://', 'https://'), 'chr21', 0, 100)
    assert np.allclose(x_local, x_remote, equal_nan=True)

    x_local = bbi.fetch(BB_FILE, 'chr21', 0, 100)
    x_remote = bbi.fetch(
        BB_URL.replace('http://', 'https://'), 'chr21', 0, 100)
    assert np.allclose(x_local, x_remote, equal_nan=True)


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_missing(path):
    x = bbi.fetch(path, 'chr21', 0, 1000, oob=0)
    assert np.all(x[:10] == 0)
    x = bbi.fetch(path, 'chr21', 0, 1000, missing=np.nan)
    assert np.all(np.isnan(x[:10]))


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_oob(path):
    x = bbi.fetch(path, 'chr21', -10, 1000, oob=np.nan)
    assert np.all(np.isnan(x[:10]))
    x = bbi.fetch(path, 'chr21', -10, 1000, oob=0)
    assert np.all(x[:10] == 0)

    n = bbi.chromsizes(path)['chr21']
    x = bbi.fetch(path, 'chr21', n - 1000, n + 10, oob=np.nan)
    assert np.all(np.isnan(x[-10:]))
    x = bbi.fetch(path, 'chr21', n - 1000, n + 10, oob=0)
    assert np.all(x[-10:] == 0)


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_intervals(path):
    x = list(bbi.fetch_intervals(path, 'chr21', 0, 1000))  # unmappable region
    assert len(x) == 0
    x = list(bbi.fetch_intervals(path, 'chr21', 0, 10000000))
    assert len(x) > 0


@pytest.mark.parametrize('path', bbi_paths)
def test_fetch_summary_stats(path):
    x = bbi.fetch(path, 'chr21', 20000000, 20001000, bins=10, summary='mean')
    y = bbi.fetch(path, 'chr21', 20000000, 20001000, bins=10)
    assert np.allclose(x, y)

    values = bbi.fetch(path, 'chr21', 20000000, 20001000)
    vmin = bbi.fetch(
        path, 'chr21', 20000000, 20001000, bins=10, summary='min'
    ).min()
    assert np.isclose(vmin, np.min(values))
    vmax = bbi.fetch(
        path, 'chr21', 20000000, 20001000, bins=10, summary='max'
    ).max()
    assert np.isclose(vmax, np.max(values))

    with pytest.raises(ValueError):
        bbi.fetch(path, 'chr21', 20000000, 20001000, bins=10, summary='foo')


def test_follows_all_redirect_codes_and_bypasses_head403():
    # Upstream kent source doesn't check for redirect status codes 307 and 308
    # So this would fail had we not patched it
    url = 'https://www.encodeproject.org/files/ENCFF620UMO/@@download/ENCFF620UMO.bigWig'
    bbi.fetch(url, 'chr21', 0, 1000)


def test_manual_redirect():
    # For some reason, HIPPAA-compliant signed AmazonS3 URLs generate 403 when
    # following redirects via HEAD (but not via GET)
    # See https://stat.ethz.ch/pipermail/bioc-devel/2016-May/009241.html
    url = 'https://www.encodeproject.org/files/ENCFF620UMO/@@download/ENCFF620UMO.bigWig'
    r = requests.head(url, allow_redirects=True)
    assert r.status_code == 403
