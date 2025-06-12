#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <span>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

/**
 * @brief Computes a convex polytope matching an extended Gaussian image.
 *
 * This function reconstructs a convex polytope whose facet normals and areas correspond to an extended Gaussian image,
 * as described in Little, James J. "An iterative method for reconstructing convex polyhedra from extended Gaussian
 * images." https://dl.acm.org/doi/proceedings/10.5555/2886844
 *
 * @param normals A span of 3D vectors representing the directions of the polytope's facets.
 * @param areas A span of doubles specifying the area of each corresponding facet.
 * @param filename A filename to which the resulting polytope is written in case the optimization is terminated
 * prematurely.
 *
 * @return A triangulated convex polytope represented as a CGAL Polyhedron.
 *
 * @note The sizes of the normals and areas spans must match.
 */
CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel>
compute_caustic(std::span<std::array<double, 3>> normals, std::span<double> areas, const char *filename);
