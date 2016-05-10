.PHONY: all build clean clean-py clean-c libkent

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
export LDFLAGS=-L/net/proteome/home/nezar/local/devel/bbifile/pykent/lib -L/usr/lib -lz -lc -lpthread
export INC=-I/net/proteome/home/nezar/local/devel/bbifile/pykent/include -I/net/proteome/home/nezar/local/devel/bbifile/pykent/src -I/usr/include
export DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}


all: build

build: lib/$(MACHTYPE)/libkent.a
	python setup.py build_ext --inplace

lib/$(MACHTYPE)/libkent.a: libkent
	mkdir -p lib/$(MACHTYPE)
	ar rcus lib/$(MACHTYPE)/libkent.a src/*.o

libkent:
	cd src && $(MAKE)

clean-c:
	rm -f lib/$(MACHTYPE)/libkent.a	
	cd src && rm -f *.o

clean-py:
	rm pykent/kent.c pykent/*.so
	rm -rf build

clean: clean-py clean-c

