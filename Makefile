
test: masktovtk pfsoltovtk
	./masktovtk mask-test-1.txt mask-test-1.vtk
	./masktovtk mask-test-2.txt mask-test-2.vtk
	./masktovtk mask-test-3.txt mask-test-3.vtk
	./masktovtk mask-test-4.txt mask-test-4.vtk

clean:
	rm -f *.vtk
	rm -f masktovtk
	rm -f pfsoltovtk

masktovtk: masktovtk.cpp simplify.h
	g++ -std=c++11 masktovtk.cpp -o masktovtk

pfsoltovtk: pfsoltovtk.c
	gcc  pfsoltovtk.c -o pfsoltovtk
