
#include <stdio.h>

int main(int argc, char** argv)
{
   int i, j;
   FILE* inFile;
   FILE* outFile;
   char* inFilename=argv[1];
   char* outFilename=argv[2];
   int version;
   int numVertices;
   int numSolids;
   int numTriangles;
   int numPatches;

   inFile = fopen(inFilename, "r");
   outFile = fopen(outFilename, "w");

   fprintf(outFile, "# vtk DataFile Version 2.0\n");
   fprintf(outFile, "%s\n", inFilename);
   fprintf(outFile, "ASCII\n");

   fscanf(inFile, "%d", &version);

   fscanf(inFile, "%d", &numVertices);

   fprintf(outFile, "DATASET POLYDATA\n");
   fprintf(outFile, "POINTS %d float\n", numVertices);

   for(i = 0; i < numVertices; ++i)
   {
      double x, y, z;
      fscanf(inFile, "%lf %lf %lf", &x, &y, &z);
      fprintf(outFile, "%f %f %f\n", x, y, z);
   }

   fscanf(inFile, "%d", &numSolids);

   for(i = 0; i < numSolids; ++i)
   {
      fscanf(inFile, "%d", &numTriangles);
      fprintf(outFile, "POLYGONS %d %d\n", numTriangles, (3+1)*numTriangles);

      for(j = 0; j < numTriangles; ++j)
      {
	 int p1, p2, p3;
	 fscanf(inFile, "%d %d %d", &p1, &p2, &p3);
	 fprintf(outFile, "%d %d %d %d\n", 3, p1, p2 , p3);
      }
   }

   fscanf(inFile, "%d", &numPatches);
   printf("NumPatches = %d\n", numPatches);

   for(i = 0; i < numPatches; ++i)
   {
      fscanf(inFile, "%d", &numTriangles);
      printf("\tNumTriangles = %d\n", numTriangles);
      for(j = 0; j < numTriangles; ++j)
      {
	 int t;
	 fscanf(inFile, "%d", &t);
      }
   }

   fclose(inFile);
}
