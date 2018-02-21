
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

using namespace std;

int main(int argc, char **argv)
{
  ifstream mask("mask-test-1.txt");

  int nx, ny, nz;

  mask >> nx >> ny >> nz;

  cout << "NX = " << nx << " NY =" << ny << " NZ = " << nz << std::endl;

  assert(nz == 1);

  vector<char> index(nx*ny);

  for(int i = 0; i < nx; ++i)
  {
    for(int j = 0; j < ny; ++j)
    {
      int ind;
      mask >> ind;
      index[ j * ny + i ] = ind;
    }
  }

  for(int i = 0; i < nx; ++i)
  {
    for(int j = 0; j < ny; ++j)
    {
    }
  }
}
