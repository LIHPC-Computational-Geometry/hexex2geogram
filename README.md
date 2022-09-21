# hexex2geogram

A `.hexex` file stores a tetrahedral mesh with an $u,v,w$ parametrization. The syntax is explained in the additional material of [HexEx](https://www.graphics.rwth-aachen.de/publication/03260/)[^HexEx] and also [below](#hexex-format).

`.geogram` is a file format used by the [Geogram library](https://github.com/BrunoLevy/geogram/wiki/Mesh), which supports embedded attributes.

This simple program converts a `.hexex` file to a `.geogram` one, using per-cell-corner attributes to store the parametrization.

## Requirements

- a C++ 17 compiler
- [CMake](https://cmake.org/)
- [UltiMaille](https://github.com/ssloy/ultimaille) (included as submodule)
- [Graphite](https://github.com/BrunoLevy/GraphiteThree) or Vorpaview, to visualize the output `.geogram` file
- [OpenMP](https://www.openmp.org) (optionnal)

## Build

```bash
mkdir build
cd build
cmake ..
make
```

## Test

The `sphere.hexex` file comes from the additional material of [HexEx](https://www.graphics.rwth-aachen.de/publication/03260/)[^HexEx].

```bash
./hexex2geogram ../data/sphere.hexex ../data/sphere.geogram
```

Then:
- open `sphere.geogram` with [Graphite](https://github.com/BrunoLevy/GraphiteThree)
- in the Properties panel, choose "ATTRIBUTE" for "painting"
- choose "cell_corners.u", "cell_corners.v" or "cell_corners.w" for "attribute"
- click on "autorange"
- use a perceptually correct colormap

The output files of [MC3D](https://github.com/HendrikBrueckler/MC3D)[^MC3D] also use the `.hexex` format, with additional information at the end (the walls of the block decomposition). These walls could be saved in the output `.geogram` file with cell facets attributes, but as of today (September 2022), [Graphite/Vorpaview cannot display them](https://github.com/BrunoLevy/geogram/issues/19). Instead, tetrahedra are grouped by block, and a cell attribute "cells.block_id" is exported. The computation is quite slow.

## `.hexex` format

It is an ASCII-based format.

The fist line is $n$, the number of vertices. It is followed by $n$ lines defining the $n$ vertices. A vertex definition has 3 floating-point numbers for the $x$, $y$ and $z$ coordinates, separated by spaces.

Then there is $m$, the number of cells, followed by $m$ lines for the tetrahedra definitions. A tetrahedron defintion has 4 integers $i0,i1,i2,i3$ (vertex index for each corner) then $4 \times 3$ floating-point numbers ( $u,v,w$ for corner 0, then for corner 1, corner 2 and corner 3), separated by spaces.

The vertices indices $i0,i1,i2,i3$ must be ordered such that
$$\det(i1-i0,i2-i0,i3-i0)>0$$
This condition is not checked in hexex2geogram.

In case of an output of [MC3D](https://github.com/HendrikBrueckler/MC3D)[^MC3D], the file ends with the number of wall triangles $w$, followed by $w$ lines for wall triangle definitions. A wall triangle definition has 3 integers and a floating-point number (separated by spaces): the 3 vertex indices and the distance from the brush fire's origin.

[^HexEx]:
    Max Lyon, David Bommes, Leif Kobbelt, *HexEx: Robust Hexahedral Mesh Extraction*, SIGGRAPH 2016, [url](https://www.graphics.rwth-aachen.de/publication/03260/)
[^MC3D]:
    Hendrik Brückler, Ojaswi Gupta, Manish Mandad, Marcel Campen, *The 3D Motorcycle Complex for Structured Volume Decomposition*, Eurographics 2022, [url](http://graphics.cs.uos.de/papers/3D_Motorcycle_Graph_EG2022.pdf)