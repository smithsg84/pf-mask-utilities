/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2018, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 **********************************************************************EHEADER*/

#include "simplify.h"
#include "error.h"
#include "readdatabox.h"

#include "tclap/CmdLine.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

/* Function IsValidFileType - This function is used to make sure a given file */
/* type is a valid one.                                                       */
/*                                                                            */
/* Parameters                                                                 */
/* ----------                                                                 */
/* char *option - The file type option to be validated                        */
/*                                                                            */
/* Return value - int - One if the file type option is valid and zero other-  */
/*                      wise                                                  */

int IsValidFileType(
                    char *option)
{
  if (strcmp(option, "pfb") == 0
      || strcmp(option, "pfsb") == 0
      || strcmp(option, "sa") == 0
      || strcmp(option, "sa2d") == 0     // Added @ IMF
      || strcmp(option, "sb") == 0
      || strcmp(option, "fld") == 0
      || strcmp(option, "vis") == 0
      || strcmp(option, "asc") == 0
#ifdef HAVE_SILO
      || strcmp(option, "silo") == 0
#endif
      || strcmp(option, "rsa") == 0)
    return(1);
  else
    return(0);
}


/* Function GetValidFileExtension - This function is used to determine the    */
/* extension of any given file name and determine if it is valid.             */
/*                                                                            */
/* Parameters                                                                 */
/* ----------                                                                 */
/* char *filename - The filename whose extension will be determined           */
/*                                                                            */
/* Return value - char * - A valid file extension or null if the file's       */
/*                         extension was invalid                              */

char *GetValidFileExtension(
char *filename)
{
  char *extension;

  /* Point the last character of the string */
  extension = filename + (strlen(filename) - 1);

  while (*extension != '.' && extension != filename)
    extension--;

  extension++;

  if (IsValidFileType(extension))
    return(extension);

  else
    return(NULL);
}

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

  // Write out patch labeling
  vtkFile << "CELL_DATA " << triangles -> size() << std::endl;
  vtkFile << "SCALARS patch_index int 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for (auto it = triangles -> begin(); it != triangles -> end(); ++it)
  {
    vtkFile << (*it).patch << std::endl;
  }

  vtkFile.close();
}


Databox *loadFile(char* filename)
{
  double default_value = 0.0;
  Databox    *databox;

  char* filetype;
  /* Make sure the file extension is valid */

  if ((filetype = GetValidFileExtension(filename)) == (char*)NULL)
  {
    std::cerr << "Invalid file extension on filename : " << filename << std::endl;
    exit(-1);
  }
  
  if (strcmp(filetype, "pfb") == 0)
    databox = ReadParflowB(filename, default_value);
  else if (strcmp(filetype, "pfsb") == 0)
    databox = ReadParflowSB(filename, default_value);
  else if (strcmp(filetype, "sa") == 0)
    databox = ReadSimpleA(filename, default_value);
  else if (strcmp(filetype, "sb") == 0)
    databox = ReadSimpleB(filename, default_value);
  else if (strcmp(filetype, "fld") == 0)
    databox = ReadAVSField(filename, default_value);
  else if (strcmp(filetype, "silo") == 0)
    databox = ReadSilo(filename, default_value);
  else if (strcmp(filetype, "asc") == 0)
    databox = ReadASCMask(filename, default_value);
  else
    databox = ReadRealSA(filename, default_value);

  return databox;

}

int main(int argc, char **argv)
{

  string inFilename;
  string vtkOutFilename;
  string pfsolOutFilename;
  int bottom;
  int side;
  float depth;

  try {  

    // Define the command line object.
    TCLAP::CmdLine cmd("Convert mask files to pfsol file", ' ', "1.0");

    TCLAP::ValueArg<string> inFilenameArg("m","mask","Mask filename",true,"mask.pfb","string");
    cmd.add( inFilenameArg );

    TCLAP::ValueArg<string> vtkOutFilenameArg("v","vtk","VTK ouput filename",true,"output.vtk","string");
    cmd.add( vtkOutFilenameArg );

    TCLAP::ValueArg<string> pfsolOutFilenameArg("s","pfsol","PFSOL ouput filename",true,"output.pfsol","string");
    cmd.add( pfsolOutFilenameArg );

    TCLAP::ValueArg<int> bottomArg("b","bottom","Bottom index",true,10,"int");
    cmd.add( bottomArg );

    TCLAP::ValueArg<int> sideArg("i","side","Side index",true,11,"int");
    cmd.add( sideArg );

    TCLAP::ValueArg<float> depthArg("d","depth","Override depth from mask file",false,NAN,"float");
    cmd.add( depthArg );

    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg. 
    inFilename = inFilenameArg.getValue();;
    vtkOutFilename = vtkOutFilenameArg.getValue();
    pfsolOutFilename = pfsolOutFilenameArg.getValue();
    bottom = bottomArg.getValue();
    side = sideArg.getValue();;
    depth = depthArg.getValue();;
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  { 
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
  }

  int nx, ny, nz;
  double sx = 0, sy = 0, sz = 0;
  double dx, dy, dz;

  Databox    *databox;

  char* c_filename = strdup(inFilename.c_str());

  // patch_names = GetStringDefault(key,
  // "left right front back bottom top");

  databox = loadFile(c_filename);

  nx = DataboxNx(databox);
  ny = DataboxNy(databox);
  nz = DataboxNz(databox);

  dx = DataboxDx(databox);
  dy = DataboxDy(databox);

  if(isnan(depth))
  {
    dz = DataboxDz(databox);
  }
  else
  {
    dz = depth;
  }
  cout << "Domain Size = (" << nx << "," << ny << "," << nz << ")" << std::endl;
  cout << "Cell Size = (" << dx << "," << dy << "," << dz << ")" << std::endl;
  cout << "Bottom patch = " << bottom << std::endl;
  cout << "Side patch = " << side << std::endl;

  assert(nz == 1);

  vector<int> indicators(nx*ny);

  for(int j = 0; j < ny; ++j)
  {
    for(int i = 0; i < nx; ++i)
    {
      double indicator;
      int k = 0;
      indicators[ triangleIndex(i,j,0) ] = *(DataboxCoeff(databox, i, j, k));
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

      if (indicator != 0 )
      {
	// Top
	{
	  Simplify::Triangle triangle;
	  triangle.patch = indicator;

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
	  triangle.patch = bottom;

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
	  triangle.patch = side;

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
	  triangle.patch = side;

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
	  triangle.patch = side;

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
	  triangle.patch = side;

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
