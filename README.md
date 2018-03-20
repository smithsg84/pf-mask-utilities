# ParFlow Mask Utilities

## Building

Requires C++ compiler.

make

## Testing

make test

## mask-to-pfsol

Utitility to build 3D PFSOL domain from 2D mask file.

### Usage

mask-to-pfsol <mask input filename> <VTK output filename> <PFSOL output filename> <bottom patch id> <side_patch id>

Creates a PFSOL file based on 2D mask input that defines the domain
and patches on the top surface.  The domain is extruded in Z to form a
3D domain.  The mask input file can be many of the standard file types
supported by ParFlow, such as ParFlow binary or simple ASCI. The ASC
file format is also supported.

The bottom and side patches are labeled with the supplied patch id's.

The mask input must be 2D with number of points in Z = 1;

The VTK file output is also generated for easier visualization.  The
VTK output includes the patch indices as data on the triangles making
the patches viewable in many tools (such as Paraview, Visit).

### ASC file format

The input mask is an ASC file format with the following format:

ncols        4
nrows        4
xllcorner    0.0
yllcorner    0.0
cellsize     1.0
NODATA_value  0.0
<ncols * nrows values>

A 0 value is outside the domain, any other value is inside the domain.

## pfsol-to-vtk

This utility is used to convert a PFSOL file to a VTK for easier visualization.

### Usage

pfsol-to-vtk <PFSOL input filename> <VTK output filename>




