.PHONY: all build clean build-cython clean-cython build-ucsc clean-ucsc test build-sdist publish-test publish

current_dir = $(shell pwd)
UNAME_S := $(shell uname -s)

ifeq (${MACHTYPE},)
	MACHTYPE:=$(shell uname -m)
#    $(info MACHTYPE was empty, set to: ${MACHTYPE})
endif
ifneq (,$(findstring -,$(MACHTYPE)))
#    $(info MACHTYPE has - sign ${MACHTYPE})
	MACHTYPE:=$(shell uname -m)
#    $(info MACHTYPE has - sign set to: ${MACHTYPE})
endif
export MACHTYPE

export CC=gcc
export COPTS=-g -pthread -fPIC -static
export CFLAGS=-Wall
export LDFLAGS=-L${current_dir}/src/${MACHTYPE} -L${current_dir}/htslib -L/usr/lib -lz -lc -lpthread -lssl -lcrypto -lpng -lhts
export INC=-I${current_dir}/include -I${current_dir}/htslib -I/usr/include
export DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE} -DUSE_SSL


ifeq (${SETUP_PY},)
	# find openssl
	SSLFLAGS=
	SLLINCL=
	ifeq (${SSLDIR},)
		SSLDIR=/usr/include/openssl
	endif
	ifneq (${SSLDIR}, "/usr/include/openssl")
		ifneq ($(UNAME_S),Darwin)
			SSLFLAGS=-L${SSLDIR}/lib
		endif
		SSLINCL=-I${SSLDIR}/include
	endif
	LDFLAGS+=${SSLFLAGS}
	INC+=${SSLINCL}

	# find libpng
	ifeq (${PNGLIB},)
		ifneq ($(wildcard /usr/lib64/libpng.a),)
		  PNGLIB=-L/usr/lib64/libpng.a
		endif
	endif
	ifeq (${PNGLIB},)
		ifneq ($(wildcard /usr/lib/libpng.a),)
		  PNGLIB=-L/usr/lib/libpng.a
		endif
	endif
	ifeq (${PNGLIB},)
		ifneq ($(wildcard /opt/local/lib/libpng.a),)
		  PNGLIB=-L/opt/local/lib/libpng.a
		endif
	endif
	ifeq (${PNGLIB},)
		ifneq ($(wildcard /usr/local/lib/libpng.a),)
		  PNGLIB=-L/usr/local/lib/libpng.a
		endif
	endif
	ifeq (${PNGLIB},)
		PNGLIB := $(shell libpng-config --ldflags  || true)
	endif
	ifeq (${PNGLIB},)
		PNGLIB=-lpng
	endif
	ifeq (${PNGINCL},)
		ifneq ($(wildcard /opt/local/include/png.h),)
		  PNGINCL=-I/opt/local/include
		else
		  PNGINCL := $(shell libpng-config --I_opts  || true)
		endif
	endif
	LDFLAGS+=${PNGLIB}
	INC+=${PNGINCL}
endif

# pass through COREDUMP
ifneq (${COREDUMP},)
	DEFS+=-DCOREDUMP
endif


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


test:
	pytest

build-sdist: clean
	python setup.py sdist

# pip install --index-url https://test.pypi.org/simple/ pybbi
publish-test: build-sdist
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

publish: build-sdist
	twine upload dist/*
