# ParFlow Mask Utilities

### Building

make

### Testing

make test

### Usage

mask-to-pfsol <mask input filename> <VTK output filename> <PFSOL output filename>

The mask input can be the standard file types supported by ParFlow, such as ParFlow binary, simple asci.

The mask input must be 2D with number of points in Z = 1;

The ASC file format is also supported.

The VTK file format includes the patch indices as data on the
triangles making the patches viewable in many tools (such as Paraview, Visit).

## ASC file format

The input mask is an ASC file format with the following format:

ncols        4
nrows        4
xllcorner    0.0
yllcorner    0.0
cellsize     1.0
NODATA_value  0.0
<ncols * nrows values>

A 0 value is outside the domain, any other value is inside the domain.






