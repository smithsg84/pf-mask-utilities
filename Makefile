
test: masktovtk pfsoltovtk maskdownsize
	./masktovtk mask-test-1.txt mask-test-1.vtk
	./masktovtk mask-test-2.txt mask-test-2.vtk
	./masktovtk mask-test-3.txt mask-test-3.vtk
	./masktovtk mask-test-4.txt mask-test-4.vtk
	./masktovtk mask-test-5.txt mask-test-5.vtk
	./masktovtk mask-test-6.txt mask-test-6.vtk

test2: masktovtk pfsoltovtk maskdownsize
	./maskdownsize mask.NA_1km.txt test.txt
	./masktovtk test.txt test.vtk
	paraview test.vtk 

clean:
	rm -f *.vtk
	rm -f masktovtk
	rm -f pfsoltovtk
	rm -f maskdownsize

masktovtk: masktovtk.cpp simplify.h
	g++ -O3 -std=c++11 masktovtk.cpp -o masktovtk

maskdownsize: maskdownsize.cpp
	g++ -O3 -std=c++11 maskdownsize.cpp -o maskdownsize

pfsoltovtk: pfsoltovtk.c
	gcc -O3  pfsoltovtk.c -o pfsoltovtk
