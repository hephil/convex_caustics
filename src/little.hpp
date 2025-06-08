#pragma once

#include <vector>
#include <assert.h>
#include <float.h>
#include <math.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>


class ErrorFunction
{
    std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Vector_3> n;

public:
    std::vector<double> A_G;

    ErrorFunction(std::vector<std::vector<double>>, std::vector<double>);

    std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3> halfspaces(std::vector<double>);
    double f(std::vector<double> x);
    std::vector<double> G(CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel>);
};

void print_vector(std::vector<double>);
CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> construct_polytope(
    std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3>
);
void test_polytope();
CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 compute_centroid(
    CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel>
);
double compute_volume(
    CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel>,
    CGAL::Exact_predicates_inexact_constructions_kernel::Point_3
);
double vector_length(std::vector<double> &);
void scale_vector(std::vector<double> &, double);
int little(std::vector<double> &, ErrorFunction);
void test_little();
