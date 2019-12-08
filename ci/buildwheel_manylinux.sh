#!/bin/bash
# This bash script is based on https://github.com/joerick/cibuildwheel
set -o errexit
set -o xtrace

mkdir -p /tmp/located
mkdir -p /tmp/delocated
mkdir -p /tmp/wheelhouse

# install system dependencies
yum install -y gcc make zlib-devel openssl-devel libpng-devel

# install a patched patchelf for auditwheel
git clone https://github.com/nvictus/patchelf.git /patchelf
cd /patchelf && ./bootstrap.sh && ./configure && make && make install

# let's go
PY_VERSIONS='cp36-cp36m cp37-cp37m cp38-cp38'

cd /project

for PY_VERSION in $PY_VERSIONS; do
    PYTHON=/opt/python/${PY_VERSION}/bin/python
    PIP=/opt/python/${PY_VERSION}/bin/pip
    rm -rf /tmp/located* /tmp/delocated/*

    # Build the package in-place
    $PIP install numpy cython
    $PIP install -v -e .

    # Quick test
    $PYTHON -c "import bbi"
    
    # Build the wheel
    $PIP wheel -v -w /tmp/located --no-deps .

    # Repair the wheel
    built_wheel=(/tmp/located/*.whl)
    auditwheel repair -w /tmp/delocated $built_wheel
    
    # Move it to output
    repaired_wheels=(/tmp/delocated/*.whl)
    mv ${repaired_wheels[@]} /tmp/wheelhouse
done
