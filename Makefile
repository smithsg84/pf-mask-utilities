SHELL := /bin/bash

MASK_TEST_FILES=$(wildcard tests/mask-test*.asc)
MASK_TESTS=$(notdir $(sort $(subst .asc,.pfsol,$(MASK_TEST_FILES))))

DEPTH_TEST_FILES=$(wildcard tests/depth-test*.asc)
DEPTH_TESTS=$(notdir $(sort $(subst .asc,.pfsol,$(DEPTH_TEST_FILES))))

all : mask-to-pfsol pfsol-to-vtk maskdownsize

test: mask-to-pfsol pfsol-to-vtk maskdownsize
	rm -f $(MASK_TESTS)
	make $(MASK_TESTS)
	rm -f $(DEPTH_TESTS)
	make TEST_ARGS="--depth 1000.0" $(DEPTH_TESTS)

%.pfsol: tests/%.asc mask-to-pfsol
	./mask-to-pfsol --mask $< --vtk $(patsubst %.pfsol,%.vtk,$@) --pfsol $@ --bottom 2 --side 3 $(TEST_ARGS)
	cmp $@ regression-test/$@
	cmp $(patsubst %.pfsol,%.vtk,$@) regression-test/$(patsubst %.pfsol,%.vtk,$@)

clean:
	rm -f *.vtk *.pfsol
	rm -f mask-to-pfsol
	rm -f pfsol-to-vtk
	rm -f maskdownsize

SRC=mask-to-pfsol.cpp readdatabox.c databox.c tools_io.c
HEADERS=databox.h parflow_config.h readdatabox.h simplify.h tools_io.h

mask-to-pfsol: $(SRC) $(HEADERS)
	g++ -Ithird-party -O3 -Wno-write-strings -std=c++11 $(SRC) -o mask-to-pfsol

maskdownsize: maskdownsize.cpp 
	g++ -O3 -std=c++11 maskdownsize.cpp -o maskdownsize

pfsol-to-vtk: pfsol-to-vtk.c
	gcc -O3  pfsol-to-vtk.c -o pfsol-to-vtk
