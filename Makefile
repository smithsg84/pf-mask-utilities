SHELL := /bin/bash

TEST_FILES=$(wildcard tests/*.asc)
TESTS=$(notdir $(sort $(subst .asc,.pfsol,$(TEST_FILES))))

all : mask-to-pfsol pfsol-to-vtk maskdownsize

test: mask-to-pfsol pfsol-to-vtk maskdownsize
	rm -f $(TESTS)
	make $(TESTS)

%.pfsol: tests/%.asc mask-to-pfsol
	./mask-to-pfsol $< $(patsubst %.pfsol,%.vtk,$@) $@ 2 3
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
	g++ -O3 -Wno-write-strings -std=c++11 $(SRC) -o mask-to-pfsol

maskdownsize: maskdownsize.cpp 
	g++ -O3 -std=c++11 maskdownsize.cpp -o maskdownsize

pfsol-to-vtk: pfsol-to-vtk.c
	gcc -O3  pfsol-to-vtk.c -o pfsol-to-vtk
