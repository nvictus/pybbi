.PHONY: all build clean clean-py clean-c

current_dir = $(shell pwd)

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
export LDFLAGS=-L${current_dir}/src/${MACHTYPE} -L/usr/lib -lz -lc -lpthread
export INC=-I${current_dir}/include -I${current_dir}/src -I/usr/include
export DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}


all: build

build: src/$(MACHTYPE)/libkent.a
	python setup.py build_ext --inplace

build-c: src/$(MACHTYPE)/libkent.a

src/$(MACHTYPE)/libkent.a:
	cd src && $(MAKE)

clean-c:
	cd src && ${MAKE} clean

clean-py:
	rm bbi/*.c bbi/*.so
	rm -rf build

clean: clean-py clean-c

