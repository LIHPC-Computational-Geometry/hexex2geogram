# hexex2geogram

A `.hexex` file stores a tetrahedral mesh with an $u,v,w$ parametrization. The syntax is explained in the additional material of [HexEx](https://www.graphics.rwth-aachen.de/publication/03260/)[^HexEx] and also [below](#hexex-format).

`.geogram` is a file format used by the [Geogram library](https://github.com/BrunoLevy/geogram/wiki/Mesh), which supports embedded attributes.

This simple program converts an `.hexex` file to a `.geogram` one, using per-cell-corner attributes to store the parametrization.

## Requirements

- a C++ 17 compiler
- [CMake](https://cmake.org/)
- [UltiMaille](https://github.com/ssloy/ultimaille) (included as submodule)
- [Graphite](https://github.com/BrunoLevy/GraphiteThree) or Vorpaview, to visualize the output `.geogram` file

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

## `.hexex` format

It is an ASCII-based format.

The fist line is $n$, the number of vertices. It is followed by $n$ lines defining the $n$ vertices. A vertex definition has 3 floating-point numbers for the $x$, $y$ and $z$ coordinates, separated by spaces.

Then there is $m$, the number of cells, followed by $m$ lines for the tetrahedra definitions. A tetrahedron defintion has 4 integers $i0,i1,i2,i3$ (vertex index for each corner) then $4 \times 3$ floating-point numbers ($u,v,w$ for corner 0, then for corner 1, corner 2 and corner 3), separated by spaces.

The vertices indices $i0,i1,i2,i3$ must be ordered such that
$$\det(i1-i0,i2-i0,i3-i0)>0$$
This condition is not checked in hexex2geogram.

## References

[^HexEx]:
    Max Lyon, David Bommes, Leif Kobbelt, *HexEx: Robust Hexahedral Mesh Extraction*, SIGGRAPH 2016
