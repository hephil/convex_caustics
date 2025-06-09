# Convex Caustics from greyscale images

This program takes a greyscale picture and computes a convex lense which projects the input image when lit by a point light source from behind.


## Explanation

The program performs three steps to compute a convec lense from a greyscale picture:

1. Derive a target distribution of surface normals (extended gaussian image) from the greyscale image.
2. Optimise the shape of a convex Polytope to match the target normal distribution using [4].
3. Transform the convex Polytope from a half-space representation (H-description) into a vertex representation (V-description) and output the geometry to a file.

### Limitations

The method makes some simplifying assumptions and approximations. Among others:
- The light source is assumed to be directional. This leads to each pixel of the image appearing as a dot in the projected image, rather than filling out the image when the lens is lit by a small light source like a phone flash.
- A simplified refraction calculation is used to compute the surface distribution that determines the shape. This leads to some warping of the final image.
- Fresnel reflectivity and 2nd order light interactions are ignored. This leads to some minor stray lights.

## Build Instructions

Requirements:
- CMake
- Gcc or Clang
- CGAL 4.10 or higher (tested up to version 6.0.1)

Install the requirements via vcpkg (instructions here: https://learn.microsoft.com/de-de/vcpkg/get_started/overview) via:

    vcpkg install cgal

then simply build the repository with

    mkdir build/
    cd build
    cmake ..
    cmake --build .


## Usage

Once built, the resulting binary can be used like this:

./convcaust input.pgm output.off


## Examples

The examples folder contains the result of an optimisation. The file `examples/olympic_rings/lens.stl` contains the mesh that was used to CNC mill the final lens out of acrylic. For fabrication, the model was inverted to create a negative imprint. Apart from that, only minor manual modifications were made to the geometry.

### Input
![photo](examples/olympic_rings/input.png)

### Mesh

![photo](examples/olympic_rings/mesh.jpg)

### Simulation

![photo](examples/olympic_rings/simulation.jpg)

### Real Photo

![photo](examples/olympic_rings/photo.jpg)


## References

- [1] Piovarči, Michal, et al. "Directional screens." Proceedings of the 1st Annual ACM Symposium on Computational Fabrication. 2017.
- [2] Schwartzburg, Yuliy, et al. "High-contrast computational caustic design." ACM Transactions on Graphics (TOG) 33.4 (2014): 1-11.
- [3] Weyrich, Tim, et al. "Fabricating microgeometry for custom surface reflectance." ACM Transactions on Graphics (TOG) 28.3 (2009): 1-6.
- [4] Little, James J. "An iterative method for reconstructing convex polyhedra from extended Gaussian images." Proceedings of the third AAAI conference on artificial intelligence. 1983.

## Note

This is a partial rewrite of a University project done in the summer semester 2017 for the lecture "Computational Fabrication".
