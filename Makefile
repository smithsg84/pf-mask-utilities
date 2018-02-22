
test: masktopfsol pfsoltovtk maskdownsize
	./masktopfsol mask-test-1.txt mask-test-1.vtk mask-test-1.pfsol
	./masktopfsol mask-test-2.txt mask-test-2.vtk mask-test-2.pfsol
	./masktopfsol mask-test-3.txt mask-test-3.vtk mask-test-3.pfsol
	./masktopfsol mask-test-4.txt mask-test-4.vtk mask-test-4.pfsol
	./masktopfsol mask-test-5.txt mask-test-5.vtk mask-test-5.pfsol
	./masktopfsol mask-test-6.txt mask-test-6.vtk mask-test-6.pfsol

test2: masktopfsol pfsoltovtk maskdownsize
	./maskdownsize mask.NA_1km.txt test.txt
	./masktopfsol test.txt test.vtk
	paraview test.vtk 

clean:
	rm -f *.vtk *.pfsol
	rm -f masktopfsol
	rm -f pfsoltovtk
	rm -f maskdownsize

masktopfsol: masktopfsol.cpp simplify.h
	g++ -O3 -std=c++11 masktopfsol.cpp -o masktopfsol

maskdownsize: maskdownsize.cpp
	g++ -O3 -std=c++11 maskdownsize.cpp -o maskdownsize

pfsoltovtk: pfsoltovtk.c
	gcc -O3  pfsoltovtk.c -o pfsoltovtk
