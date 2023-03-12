.PHONY: all build clean build-cython clean-cython build-ucsc clean-ucsc sdist publish-test publish

current_dir = $(shell pwd)
UNAME_S := $(shell uname -s)

ifeq (${MACHTYPE},)
	MACHTYPE:=$(shell uname -m)
endif
export MACHTYPE

export CC ?= gcc
export COPTS=-g -pthread -fPIC -static
export CFLAGS=-Wall $(shell pkg-config --static --cflags-only-other openssl zlib libpng)
export DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE} -DUSE_SSL
export LDFLAGS=-L${current_dir}/src/${MACHTYPE} $(shell pkg-config --static --libs openssl zlib libpng)
export INC=-I${current_dir}/include $(shell pkg-config --static --cflags-only-I openssl zlib libpng)
# pass through COREDUMP
ifneq (${COREDUMP},)
	DEFS+=-DCOREDUMP
endif

# Append values to CFLAGS and LDFLAGS
CFLAGS += $(shell echo $$CFLAGS)
LDFLAGS += $(shell echo $$LDFLAGS)


all: build

src/$(MACHTYPE)/libkent.a:
	cd src && $(MAKE)


clean-ucsc:
	cd src && ${MAKE} clean

build-ucsc: src/$(MACHTYPE)/libkent.a


clean-cython:
	rm -f bbi/*.c bbi/*.so
	find . -name '*.pyc' -exec rm --f {} +
	find . -name '*.pyo' -exec rm --f {} +

build-cython: src/$(MACHTYPE)/libkent.a
	python setup.py build_ext --inplace


clean: clean-ucsc clean-cython
	rm -rf build/
	rm -rf dist/

build: build-ucsc build-cython


sdist: clean
	python setup.py sdist

# pip install --index-url https://test.pypi.org/simple/ pybbi
publish-test: sdist
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

publish: sdist
	twine upload dist/*
