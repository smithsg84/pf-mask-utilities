
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
  string outFilename(argv[2]);

  ifstream mask(inFilename);

  int nx, ny, nz;
  double sx = 0, sy = 0, sz = 0;
  //double dx = 1000.0, dy = 1000.0 , dz = 1000.0;
  double dx = 1.0, dy = 1.0 , dz = 1.0;

  mask >> nx >> ny >> nz;

  cout << "NX = " << nx << " NY =" << ny << " NZ = " << nz << std::endl;

  assert(nz == 1);

  vector<char> indicators(nx*ny);

  for(int j = 0; j < ny; ++j)
  {
    for(int i = 0; i < nx; ++i)
    {
      int indicator;
      mask >> indicator;
      indicators[ triangleIndex(i,j,0) ] = indicator;
    }
  }

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

      if (indicator == 2 )
      {

#if 1
	// Top
	{
	  Simplify::Triangle triangle;
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
#endif

#if 0
	// Bottom
	{
	  Simplify::Triangle triangle;
	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i+1,j,0);
	  triangle.v[2]=  vertexIndex(i+1,j+1,0);
	  
	  triangles -> push_back(triangle);
	  
	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i+1,j+1,0);
	  triangle.v[2]=  vertexIndex(i,j+1,0);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j,0)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,0)].used = true;
	}

	// Left
	if ( (i == 0) || (indicators[ triangleIndex(i-1,j,0) ] == 1) )
	{
	  Simplify::Triangle triangle;
	  triangle.v[0] = vertexIndex(i,j,0);
	  triangle.v[1]=  vertexIndex(i,j,1);
	  triangle.v[2]=  vertexIndex(i,j+1,0);
	  
	  triangles -> push_back(triangle);

	  triangle.v[0] = vertexIndex(i,j,1);
	  triangle.v[1]=  vertexIndex(i,j+1,0);
	  triangle.v[2]=  vertexIndex(i,j+1,1);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j,0)].used = true;
	  (*vertices)[ vertexIndex(i,j,1)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,1)].used = true;
	}

	// Right
	if ( (i == (nx - 1)) || (indicators[ triangleIndex(i+1,j,0) ] == 1) )
	{
	  Simplify::Triangle triangle;
	  triangle.v[0] = vertexIndex(i+1,j,0);
	  triangle.v[1]=  vertexIndex(i+1,j,1);
	  triangle.v[2]=  vertexIndex(i+1,j+1,0);
	  
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
	if ( (j==0) || (indicators[ triangleIndex(i,j-1,0) ] == 1) )
	{
	  Simplify::Triangle triangle;
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
	if ( (j == (ny - 1)) || (indicators[ triangleIndex(i,j+1,0) ] == 1) )
	{
	  Simplify::Triangle triangle;
	  triangle.v[0] = vertexIndex(i,j+1,0);
	  triangle.v[1]=  vertexIndex(i+1,j+1,0);
	  triangle.v[2]=  vertexIndex(i+1,j+1,1);
	  
	  triangles -> push_back(triangle);

	  triangle.v[0] = vertexIndex(i,j+1,0);
	  triangle.v[1]=  vertexIndex(i+1,j+1,1);
	  triangle.v[2]=  vertexIndex(i,j+1,1);

	  triangles -> push_back(triangle);

	  (*vertices)[ vertexIndex(i,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i,j+1,1)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,0)].used = true;
	  (*vertices)[ vertexIndex(i+1,j+1,1)].used = true;
	}
#endif
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

  //writeVTK("test-1.vtk", new_vertices, triangles);

  Simplify::swap(*new_vertices, *triangles);

  //Simplify::simplify_mesh_lossless();

  writeVTK(outFilename, &Simplify::vertices, &Simplify::triangles);
}
