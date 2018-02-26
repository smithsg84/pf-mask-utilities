
# ParFlow Mask Utilities

## ascmask-to-pfsol

This utility is used to convert a 2D mask file describing a domain into a PFSOL file format.

Currently only 3 patches are created, top, bottom and sides.

The input mask is an ASC file format with the following format:

ncols        4
nrows        4
xllcorner    0.0
yllcorner    0.0
cellsize     1.0
NODATA_value  0.0
<ncols * nrows values>

A 0 value is outside the domain, a 1 value is inside the domain.

### Building

make

### Testing

make test

### Usage

ascmask-to-pfsol <asc mask input filename> <VTK output filename> <PFSOL output filename>

## pfsol-to-vtk

This utility is used to convert a PFSOL file to a VTK for easier visualization.

### Building

make

### Usage

pfsol-to-vtk <PFSOL input filename> <VTK output filename>





