SHELL := /bin/bash

MASK_TEST_FILES=$(wildcard tests/mask-test*.asc)
MASK_TESTS=$(notdir $(sort $(subst .asc,.pfsol,$(MASK_TEST_FILES))))

DEPTH_TEST_FILES=$(wildcard tests/depth-test*.asc)
DEPTH_TESTS=$(notdir $(sort $(subst .asc,.pfsol,$(DEPTH_TEST_FILES))))

MULTIPLE_MASK_TESTS=tests/multiple-mask-test-1

all : mask-to-pfsol pfsol-to-vtk maskdownsize

test: mask-to-pfsol pfsol-to-vtk maskdownsize
	rm -f $(MASK_TESTS)
	make $(MASK_TESTS)
	rm -f $(DEPTH_TESTS)
	make TEST_ARGS="--depth 1000.0" $(DEPTH_TESTS)
	rm -f  multiple-mask-test-1.vtk multiple-mask-test-1.pfsol
	make multiple-masks-tests

%.pfsol: tests/%.asc mask-to-pfsol
	./mask-to-pfsol --mask $< --vtk $(patsubst %.pfsol,%.vtk,$@) --pfsol $@ --bottom-patch-label 2 --side-patch-label 3 $(TEST_ARGS)
	cmp $@ regression-test/$@
	cmp $(patsubst %.pfsol,%.vtk,$@) regression-test/$(patsubst %.pfsol,%.vtk,$@)

multiple-masks-tests: 
	./mask-to-pfsol \
	--mask-top tests/multiple-mask-test-1-top.asc \
	--mask-bottom tests/multiple-mask-test-1-bottom.asc \
	--mask-left tests/multiple-mask-test-1-left.asc \
	--mask-right tests/multiple-mask-test-1-right.asc \
	--mask-front tests/multiple-mask-test-1-front.asc \
	--mask-back tests/multiple-mask-test-1-back.asc \
	--vtk multiple-mask-test-1.vtk \
	--pfsol multiple-mask-test-1.pfsol
	cmp multiple-mask-test-1.pfsol regression-test/multiple-mask-test-1.pfsol
	cmp multiple-mask-test-1.vtk regression-test/multiple-mask-test-1.vtk

clean:
	rm -f *.vtk *.pfsol
	rm -f mask-to-pfsol
	rm -f pfsol-to-vtk
	rm -f maskdownsize

SRC=mask-to-pfsol.cpp readdatabox.c databox.c tools_io.c
HEADERS=databox.h parflow_config.h readdatabox.h simplify.h tools_io.h

mask-to-pfsol: $(SRC) $(HEADERS)
	g++ -Ithird-party -g -Wno-write-strings -std=c++11 $(SRC) -o mask-to-pfsol

maskdownsize: maskdownsize.cpp 
	g++ -O3 -std=c++11 maskdownsize.cpp -o maskdownsize

pfsol-to-vtk: pfsol-to-vtk.c
	gcc -O3  pfsol-to-vtk.c -o pfsol-to-vtk
