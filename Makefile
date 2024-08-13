# Compiler and complilation flag are included in Makefile.inc
include Makefile.inc

lib:
	mkdir -p include;mkdir -p lib
	cd src; mkdir -p obj; make all

example:
	cd examples; make all

all:
	mkdir -p include; mkdir -p lib
	cd src; mkdir -p obj; make all

clean:
	cd src; make clean