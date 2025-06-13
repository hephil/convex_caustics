#include "common.hpp"
#include "little.hpp"

#include <CGAL/Cartesian.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Plane = Kernel::Plane_3;
using Direction = Kernel::Direction_3;
using Triangle = Kernel::Triangle_3;
using Tetrahedron = Kernel::Tetrahedron_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using Vertex = CGAL::Polyhedron_3<Kernel>::Vertex;
using Facet = Polyhedron::Facet;
using Halfedge_facet_circulator = Polyhedron::Halfedge_around_facet_circulator;

static const double EPSILON = 1e-9;
static const Point O = { 0.0, 0.0, 0.0 }; // origin

Polyhedron compute_caustic(std::span<std::array<double, 3>> normals, std::span<double> areas, const char *filename)
{
    Assert(normals.size() == areas.size(), "Invalid Extended Gaussian Image");

#ifdef LOGGING_ENABLED
    std::ofstream logfile(std::string(filename) + ".log");
    if (logfile.is_open())
        logfile << "Iteration,Error" << std::endl;
#endif

    // 1) initialize halfspaces with distance 1 from the origin.
    std::vector<double> L(normals.size(), 1.0);
    std::vector<double> G_L(normals.size(), 0.0);
    std::vector<double> best_L(normals.size(), 1.0);
    std::vector<double> best_G_L(normals.size(), 0.0);
    std::vector<Plane> h_desc(L.size());
    Polyhedron p;
    std::vector<Point> points;

    // normalize areas
    {
        double len2 = 0.0;
        for (size_t i = 0; i < areas.size(); i++) len2 += areas[i] * areas[i];
        len2 = std::max(EPSILON, std::sqrt(len2));
        for (size_t i = 0; i < areas.size(); i++) areas[i] /= len2;
    }

    unsigned long k = 0;
    double previous_error = DBL_MAX;
    double best_error = DBL_MAX;
    double current_error = DBL_MAX;
    double gamma = 1.0;

    while (true)
    {
        // 2) Construct P(L)
        // compute the H-representation of the convex Polytope
        for (size_t i = 0; i < L.size(); i++)
        {
            h_desc[i] = Plane(
                Point(L[i] * normals[i][0], L[i] * normals[i][1], L[i] * normals[i][2]),
                Vector{ normals[i][0], normals[i][1], normals[i][2] }
            );
        }

        // 2.1) Transform the n planes given by L into M, a set of n points in R3, using the dual transform
        p.clear();
        CGAL::halfspace_intersection_with_constructions_3(h_desc.begin(), h_desc.end(), p, Point(0, 0, 0));

        points.clear();
        for (auto iter = p.vertices_begin(); iter != p.vertices_end(); iter++) points.push_back(iter->point());

        // 2.2) Compute the convex hull of M, call it CH(M).
        // 2.3) Determine the adjacency relations of P(L) from CH(M).Calculate the locations of the vertices of P(L).
        p.clear();
        CGAL::convex_hull_3(points.begin(), points.end(), p);

        // 3) Compute the centroid of P(L).
        // Compute vertex centroid first to get valid point inside polytope
        Point vertex_centroid = { 0.0, 0.0, 0.0 };
        for (const auto &v : p.vertex_handles())
            vertex_centroid = vertex_centroid + (v->point() - O) / p.size_of_vertices();

        // Compute volume and volume centroid
        double volume = 0.0;
        Point volume_centroid = { 0.0, 0.0, 0.0 };

        for (const auto &f : p.facet_handles())
        {
            assert(f->is_triangle());

            auto h = f->halfedge();
            Point a = h->vertex()->point();
            Point b = h->next()->vertex()->point();
            Point c = h->next()->next()->vertex()->point();

            Tetrahedron t(vertex_centroid, a, b, c);
            double pyramid_volume = std::abs(t.volume());

            // Compute centroid of the tetrahedron: 0.25 * (a + b + c + d).
            Point pyramid_center = O + ((vertex_centroid - O) + (a - O) + (b - O) + (c - O)) / 4.0;

            volume += pyramid_volume;
            volume_centroid = volume_centroid + pyramid_volume * (pyramid_center - O);
        }

        volume_centroid = O + (volume_centroid - O) / volume;

        // 3) Translate the centroid of P(L) to the origin.
        Vector shift = O - volume_centroid;
        for (auto &v : p.vertex_handles()) v->point() = v->point() + shift;

        // 3) Scale L by V(L)^(1/3) to make its volume unity.
        // XXX using 1.0/3.0 leads to rounding errors which terminate the optimization early
        const double almost_one_third = 0.333;
        double scale = 1.0 / pow(volume, almost_one_third);
        for (size_t i = 0; i < L.size(); i++) L[i] *= scale;

        // scale the polytope as well so that G_L is computed correctly
        for (auto v = p.vertices_begin(); v != p.vertices_end(); v++)
        {
            Vector v2 = Vector(Point(0, 0, 0), v->point()) * scale;
            v->point() = Point(v2.x(), v2.y(), v2.z());
        }

        // 3) and the gradient of V, G(L).
        for (unsigned i = 0; i < normals.size(); i++) G_L[i] = 0.0;

        for (const auto &f : p.facet_handles())
        {
            int j = -1;
            double best = 2.0;
            Facet::Halfedge_handle h = f->halfedge();
            Point p1 = h->vertex()->point();
            Point p2 = h->next()->vertex()->point();
            Point p3 = h->next()->next()->vertex()->point();

            Triangle tri = Triangle(p1, p2, p3);
            double area = sqrt(tri.squared_area());
            if (area < EPSILON) // skip tiny facets
                continue;

            Vector v = CGAL::unit_normal(p1, p2, p3);

            for (int i = 0; i < normals.size(); i++)
            {
                double vdotn = std::abs(v * Vector(normals[i][0], normals[i][1], normals[i][2]));
                if ((1.0 - vdotn) < best)
                {
                    best = 1.0 - vdotn;
                    j = i;
                }
            }

            assert(j != -1);
            assert(best < 0.05);
            G_L[j] += area;
        }

        // normalize G_L (TODO: why???)
        double len2 = 0.0;
        for (const double &g : G_L) len2 += g * g;
        double len = std::max(EPSILON, std::sqrt(len2));
        for (auto &g : G_L) g /= len;

        // 4) Evaluate f(L);
        current_error = 0.0;
        for (size_t i = 0; i < L.size(); i++) current_error += L[i] * areas[i];

        previous_error = best_error;
        // 4) if the decrease in f is less than a pre-specified value, terminate.
        if (current_error < best_error)
        {
            if (best_error - current_error < EPSILON)
                break;
            best_error = current_error;
            for (size_t i = 0; i < L.size(); i++) best_L[i] = L[i];
            for (size_t i = 0; i < G_L.size(); i++) best_G_L[i] = G_L[i];
            gamma *= 1.1;
        }
        else
        {
            best_error = current_error;
            for (size_t i = 0; i < L.size(); i++) L[i] = best_L[i];
            for (size_t i = 0; i < G_L.size(); i++) G_L[i] = best_G_L[i];
            gamma = std::max(DBL_EPSILON, 0.5 * gamma);
        }

        // 4) Otherwise, compute a step using equation (3), update L, and repeat, starting at step 2.
        double GdotA = 0.0;
        double area_error = 0.0;
        for (size_t i = 0; i < L.size(); i++)
        {
            GdotA += areas[i] * G_L[i];
            area_error += (areas[i] - G_L[i]) * (areas[i] - G_L[i]);
        }
        area_error = std::sqrt(area_error);

        // This step is a multiple of: <A,G(L)>G(L) - A , where <x,y> is the inner product
        for (unsigned i = 0; i < L.size(); i++)
        {
            double stepping = GdotA * G_L[i] - areas[i];
            L[i] = std::max(EPSILON, L[i] + (stepping * gamma));
        }

        if (k % 100 == 0)
        {
            std::ofstream out(filename);
            if (out.is_open())
                out << p;
#ifdef LOGGING_ENABLED
            if (logfile.is_open())
                logfile << std::format("{},{}", k, area_error) << std::endl;
#endif
        }
        std::cout << std::format(
            "Iteration {:>6}, Error: {:>12.4e}, Delta: {:>12.4e}",
            k,
            area_error,
            previous_error - current_error
        ) << "\t\t\r"
                  << std::flush;
        k++;
    }

    std::cout << std::endl;

    return p;
}