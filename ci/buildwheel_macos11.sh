#!/bin/bash
# This bash script is based on https://github.com/joerick/cibuildwheel
set -o errexit
set -o xtrace

mkdir -p /tmp/located
mkdir -p /tmp/delocated
mkdir -p /tmp/wheelhouse

brew install pkg-config openssl libpng

curl -L -o /tmp/get-pip.py https://bootstrap.pypa.io/get-pip.py

export C_INCLUDE_PATH="/usr/local/include/libpng:/usr/local/opt/openssl/include:$C_INCLUDE_PATH"

PY_VERSIONS=("3.10")
PKG_URLS=(
    "https://www.python.org/ftp/python/3.10.1/python-3.10.1-macos11.pkg"
)
len=${#PKG_URLS[@]}

for (( i=0; i<$len; i++ )); do
    PY_VERSION=${PY_VERSIONS[$i]}
    PKG_URL=${PKG_URLS[$i]}
    PY_BIN="/Library/Frameworks/Python.framework/Versions/$PY_VERSION/bin"

    rm -rf /tmp/located* /tmp/delocated/*
    
    # Install python
    curl -L -o /tmp/Python.pkg $PKG_URL
    sudo installer -pkg /tmp/Python.pkg -target /

    # Install pip
    PYTHON=$PY_BIN/python3
    $PYTHON /tmp/get-pip.py --no-setuptools --no-wheel
    
    # Install wheel, delocate
    PIP=$PY_BIN/pip3
    $PIP install --upgrade setuptools wheel delocate

    # Build the package in-place
    $PIP install numpy cython
    $PIP install -v -e .

    # Quick test
    $PYTHON -c "import bbi"

    # Build the wheel
    $PIP wheel -v -w /tmp/located --no-deps .

    # Repair the wheel
    built_wheel=(/tmp/located/*.whl)
    $PY_BIN/delocate-listdeps $built_wheel && $PY_BIN/delocate-wheel -w /tmp/delocated $built_wheel

    # Move it to output
    repaired_wheels=(/tmp/delocated/*.whl)
    mv ${repaired_wheels[@]} /tmp/wheelhouse
done
