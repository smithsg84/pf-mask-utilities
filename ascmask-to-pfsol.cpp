
#include "simplify.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

#define vertexIndex(i, j, k) (((nx+1) * (ny+1) * (k)) + ((j) * (nx+1)) + (i))
#define triangleIndex(i, j, k) ((nx * ny * (k)) + ((j) * nx) + (i))

typedef std::numeric_limits< double > dbl;

bool equal(double a, double b)
{
    return fabs(a - b) < DBL_EPSILON;
}

void writePFSOL(string filename, vector<Simplify::Vertex>* vertices, vector<Simplify::Triangle>* triangles)
{
  ofstream pfsolFile(filename);
  pfsolFile.precision(dbl::max_digits10);

  // Version
  pfsolFile << "1" << std::endl;

  pfsolFile << vertices -> size() << std::endl;

  for (auto it = vertices->begin(); it != vertices->end(); ++it)
  {
    double x = (*it).p.x;
    double y = (*it).p.y;
    double z = (*it).p.z;

    pfsolFile << x << " " << y << " " << z << std::endl;
  }

  // Number of solids
  pfsolFile << "1" << std::endl;

  pfsolFile << triangles -> size() << std::endl;

  for (auto it = triangles -> begin(); it != triangles -> end(); ++it)
  {
    pfsolFile << (*it).v[0]  << " " << (*it).v[1] << " " << (*it).v[2] << std::endl;
  }

  // Number of patches
  
  // This looping is not very efficient but would require some extra data structures to eliminate.
  pfsolFile << "3" << std::endl;
  for(int patch = 1; patch < 4; ++patch)
  {
    int numTriangles = 0;
    for (auto it = triangles -> begin(); it != triangles -> end(); ++it)
    {
      if((*it).patch == patch)
      {
	numTriangles++;
      }
    }
    
    pfsolFile << numTriangles << std::endl;

    std::cout << "Number of triangles in patch " << patch << " = " << numTriangles << std::endl;

    int index = 0;
    for (auto it = triangles -> begin(); it != triangles -> end(); ++it)
    {
      if((*it).patch == patch)
      {
	pfsolFile << index << std::endl;
      }

      ++index;
    }
    
  }

  pfsolFile.close();
}

void writeVTK(string filename, vector<Simplify::Vertex>* vertices, vector<Simplify::Triangle>* triangles)
{
  ofstream vtkFile(filename);

  vtkFile.precision(dbl::max_digits10);

  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << filename << std::endl;
  vtkFile << "ASCII" << std::endl;

  vtkFile << "DATASET POLYDATA" << std::endl;
  vtkFile << "POINTS " << vertices -> size() << " float" << std::endl;

  for (auto it = vertices->begin(); it != vertices->end(); ++it)
  {
    double x = (*it).p.x;
    double y = (*it).p.y;
    double z = (*it).p.z;

    vtkFile << x << " " << y << " " << z << std::endl;
  }

  vtkFile << "POLYGONS " << triangles -> size() << " " << (3+1) * triangles -> size() << std::endl;
  for (auto it = triangles -> begin(); it != triangles -> end(); ++it)
  {
    vtkFile << "3 " <<   (*it).v[0]  << " " << (*it).v[1] << " " << (*it).v[2] << std::endl;
  }

  vtkFile.close();
}

int main(int argc, char **argv)
{
  string inFilename(argv[1]);
  string vtkOutFilename(argv[2]);
  string pfsolOutFilename(argv[3]);
  double maskValue = 1;

  cout << "Reading file " << inFilename << std::endl;

  ifstream mask(inFilename);

  int nx, ny, nz;
  double sx = 0, sy = 0, sz = 0;
  double dx = 1000.0, dy = 1000.0 , dz = 1000.0;

  string text;

  mask >> text >> nx;
  mask >> text >> ny;
  nz = 1;

  mask >> text >> text;
  mask >> text >> text;
  
  mask >> text >> dx;
  dy = dz = dx;
  
  mask >> text >> text;

  cout << "Domain Size = (" << nx << "," << ny << "," << nz << ")" << std::endl;
  cout << "Cell Size = (" << dx << "," << dy << "," << dz << ")" << std::endl;

  assert(nz == 1);

  vector<char> indicators(nx*ny);

  for(int j = 0; j < ny; ++j)
  {
    for(int i = 0; i < nx; ++i)
    {
      double indicator;
      mask >> indicator;
      // ASC files are flipped around J axis from PF ordering
      int flipped_j = (ny - 1) - j;
      indicators[ triangleIndex(i,flipped_j,0) ] = equal(indicator, maskValue);
    }
  }

  cout << endl;

  vector<Simplify::Vertex>* vertices = new vector<Simplify::Vertex>((nx+1)*(ny+1)*(nz+1));

  vector<Simplify::Triangle>* triangles = new vector<Simplify::Triangle>();

  // Build list of all possible vertices
  for(int k = 0; k < nz+1; ++k)
  {
    for(int j = 0; j < ny+1; ++j)
    {
      for(int i = 0; i < nx+1; ++i)
      {
	Simplify::Vertex *vertex = &((*vertices)[ vertexIndex(i,j,k) ]);
	vertex -> p.x = sx + i * dx;
	vertex -> p.y = sy + j * dy;
	vertex -> p.z = sz + k * dz;
	vertex -> used = false;
      }
    }
  }

  // Build triangles for faces on every cell
  for(int j = 0; j < ny; ++j)
  {
    for(int i = 0; i < nx; ++i)
    {
      int indicator = indicators[ triangleIndex(i,j,0) ];

      if (indicator == 1 )
      {
	// Top
	{
	  Simplify::Triangle triangle;
	  triangle.patch = 1;

	  triangle.v[0] = vertexIndex(i,j,1);
	  triangle.v[1]=  vertexIndex(i+1,j,1);
	  triangle.v[2]=  vertexIndex(i+1,j+1,1);
	  
	  triangles -> push_back(triangle);
	  
	  triangle.v[0] = vertexIndex(i,j,1);
	  triangle.v[1]=  vertexIndex(i+1,j+1,1);
	  triangle.v[2]=  vertexIndex(i,j+1,1);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j,1)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,1)].used = true;
	}

	// Bottom
	{
	  Simplify::Triangle triangle;
	  triangle.patch = 2;

	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i+1,j+1,0);
	  triangle.v[2]=  vertexIndex(i+1,j,0);
	  
	  triangles -> push_back(triangle);

	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i,j+1,0);	  
	  triangle.v[2]=  vertexIndex(i+1,j+1,0);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j,0)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,0)].used = true;
	}

	// Left
	if ( (i == 0) || (indicators[ triangleIndex(i-1,j,0) ] == 0) )
	{
	  Simplify::Triangle triangle;
	  triangle.patch = 3;

	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i,j,1);
	  triangle.v[2]=  vertexIndex(i,j+1,0);
	  
	  triangles -> push_back(triangle);

	  triangle.v[0] = vertexIndex(i,j,1);
	  triangle.v[1]=  vertexIndex(i,j+1,1);
	  triangle.v[2]=  vertexIndex(i,j+1,0);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j,0)].used = true;
	  (*vertices)[ vertexIndex(i,j,1)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,1)].used = true;
	}

	// Right
	if ( (i == (nx - 1)) || (indicators[ triangleIndex(i+1,j,0) ] == 0) )
	{
	  Simplify::Triangle triangle;
	  triangle.patch = 3;

	  triangle.v[0]=  vertexIndex(i+1,j+1,0);
	  triangle.v[1]=  vertexIndex(i+1,j,1);
	  triangle.v[2] = vertexIndex(i+1,j,0);
	  
	  triangles -> push_back(triangle);

	  triangle.v[0] = vertexIndex(i+1,j,1);
	  triangle.v[1]=  vertexIndex(i+1,j+1,0);
	  triangle.v[2]=  vertexIndex(i+1,j+1,1);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i+1,j,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,1)].used = true;
	}

	// Front
	if ( (j==0) || (indicators[ triangleIndex(i,j-1,0) ] == 0) )
	{
	  Simplify::Triangle triangle;
	  triangle.patch = 3;

	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i+1,j,0);
	  triangle.v[2]=  vertexIndex(i+1,j,1);
	  
	  triangles -> push_back(triangle);

	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i+1,j,1);
	  triangle.v[2]=  vertexIndex(i,j,1);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j,0)].used = true;
	  (*vertices)[ vertexIndex(i,j,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j,1)].used = true;
	}

	// Back
	if ( (j == (ny - 1)) || (indicators[ triangleIndex(i,j+1,0) ] == 0) )
	{
	  Simplify::Triangle triangle;
	  triangle.patch = 3;

	  triangle.v[2] = vertexIndex(i,j+1,0);
	  triangle.v[1]=  vertexIndex(i+1,j+1,0);
	  triangle.v[0]=  vertexIndex(i+1,j+1,1);
	  
	  triangles -> push_back(triangle);

	  triangle.v[2] = vertexIndex(i,j+1,0);
	  triangle.v[1]=  vertexIndex(i+1,j+1,1);
	  triangle.v[0]=  vertexIndex(i,j+1,1);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,1)].used = true;
	}
      }
    }
  }

  // Create new vertices vector with only used vertices.
  int new_index = 0;
  vector<Simplify::Vertex>* new_vertices = new vector<Simplify::Vertex>();

  for (auto it = vertices -> begin(); it != vertices -> end(); ++it)
  {
     if ((*it).used)
     {
	(*it).new_index = new_index++;
	new_vertices -> push_back((*it));
      }
  }

  // Reindex vertices into new compressed vertices vector.
  for (auto it = triangles -> begin(); it != triangles -> end(); ++it)
  {
    (*it).v[0] = (*vertices)[(*it).v[0]].new_index;
    (*it).v[1] = (*vertices)[(*it).v[1]].new_index;
    (*it).v[2] = (*vertices)[(*it).v[2]].new_index;
  }

  vertices = 0;

  Simplify::swap(*new_vertices, *triangles);

  Simplify::simplify_mesh_lossless();

  writeVTK(vtkOutFilename, &Simplify::vertices, &Simplify::triangles);
  writePFSOL(pfsolOutFilename, &Simplify::vertices, &Simplify::triangles);
}
