#include "common.hpp"
#include "little.hpp"

static const double EPSILON = 1e-9;
static const bool VERBOSE = false;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Direction_3 Direction;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Tetrahedron_3 Tetrahedron;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet Facet;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

ErrorFunction::ErrorFunction(std::vector<std::vector<double>> xn, std::vector<double> xA_G)
{
    for (unsigned i = 0; i < xn.size(); i++)
    {
        Vector v = Vector(xn[i][0], xn[i][1], xn[i][2]);
        this->n.push_back(v);
    }
    scale_vector(xA_G, 1 / vector_length(xA_G));
    this->A_G = xA_G;
}

std::vector<Plane> ErrorFunction::halfspaces(std::vector<double> x)
{
    assert(x.size() == n.size());
    std::vector<Plane> ret;
    for (unsigned i = 0; i < x.size(); i++)
    {
        Vector v = n[i];
        v *= x[i];
        Point p = Point(v[0], v[1], v[2]);
        Plane pl = Plane(p, v);
        ret.push_back(pl);
    }
    return ret;
}

double ErrorFunction::f(std::vector<double> L)
{
    assert(L.size() == A_G.size());
    double ret = 0.0;
    for (unsigned i = 0; i < L.size(); i++)
    {
        ret += L[i] * A_G[i];
    }
    return ret;
}

// calculates the area of facets facing in directions given by this->n.
// Assumes P is triangulated
std::vector<double> ErrorFunction::G(Polyhedron P)
{
    if (VERBOSE)
        std::cout << "Computing G(L)...";

    // initialize all with 0.0
    std::vector<double> ret;
    ret.resize(this->n.size());
    for (unsigned i = 0; i < this->n.size(); i++)
    {
        ret[i] = 0.0;
    }

    for (Polyhedron::Facet_iterator iter = P.facets_begin(); iter != P.facets_end(); iter++)
    {
        int j = -1;
        double best = 2.0;
        Facet f = *iter;
        Facet::Halfedge_handle h = f.halfedge();
        Point p1 = h->vertex()->point();
        Point p2 = h->next()->vertex()->point();
        Point p3 = h->next()->next()->vertex()->point();

        Triangle tri = Triangle(p1, p2, p3);
        double area = sqrt(tri.squared_area());
        if (area < EPSILON) // skip tiny facets
            continue;

        Vector v = CGAL::unit_normal(p1, p2, p3);

        for (int i = 0; i < this->n.size(); i++)
        {
            // Vector u = CGAL::cross_product(v, n[i]);
            if (abs((v * n[i]) - 1.0) < best)
            {
                best = abs((v * n[i]) - 1.0);
                j = i;
            }
        }

        assert(j != -1);
        assert(best < 0.05);
        ret[j] += area;
    }

    if (VERBOSE)
        std::cout << "DONE.\n";

    return ret;
}

void print_vector(std::vector<double> a)
{
    std::cout << '(';
    for (unsigned i = 0; i < a.size(); i++)
    {
        printf("%'.3f", a[i]);
        if (i + 1 != a.size())
            std::cout << ", ";
    }
    std::cout << ")";
}

Polyhedron construct_polytope(std::vector<Plane> halfspaces)
{
    if (VERBOSE)
        std::cout << "> Constructing Polytope from halfspaces...";

    std::vector<Plane>::iterator hbegin = halfspaces.begin();
    std::vector<Plane>::iterator hend = halfspaces.end();

    Polyhedron pm;
    CGAL::halfspace_intersection_with_constructions_3(hbegin, hend, pm, Point(0, 0, 0));
    Polyhedron p;

    std::vector<Point> points;
    for (auto iter = pm.vertices_begin(); iter != pm.vertices_end(); iter++) points.push_back(iter->point());

    // again convex hull but this time with unique point cloud
    CGAL::convex_hull_3(points.begin(), points.end(), p);

    if (VERBOSE)
    {
        std::cout << "DONE.\n";
    }

    return p;
}

Point compute_centroid(Polyhedron P)
{
    if (VERBOSE)
        std::cout << "> Computing centroid...";

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for (Polyhedron::Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); iter++)
    {
        Point p = iter->point();
        x += p.x();
        y += p.y();
        z += p.z();
    }
    x = x / P.size_of_vertices();
    y = y / P.size_of_vertices();
    z = z / P.size_of_vertices();
    Point geoCenter = Point(x, y, z);
    double weight = 0;
    double tmpX = 0;
    double tmpY = 0;
    double tmpZ = 0;
    for (Polyhedron::Facet_iterator iter = P.facets_begin(); iter != P.facets_end(); iter++)
    {
        Facet f = *iter;
        assert(f.is_triangle());
        auto edges = iter->halfedge();
        auto a = edges->vertex()->point();
        edges = edges->next();
        auto b = edges->vertex()->point();
        edges = edges->next();
        auto c = edges->vertex()->point();
        Tetrahedron t = Tetrahedron(geoCenter, a, b, c);
        Point bottomCenter = centroid(a, b, c);
        Point tetraCenter = centroid(t);
        Point pyramidCenter = barycenter(bottomCenter, 0.75, tetraCenter, 0.25);
        weight += abs(t.volume());
        tmpX += pyramidCenter.x() * t.volume();
        tmpY += pyramidCenter.y() * t.volume();
        tmpZ += pyramidCenter.z() * t.volume();
    }
    tmpX = tmpX / weight;
    tmpY = tmpY / weight;
    tmpZ = tmpZ / weight;
    Point center = Point(tmpX, tmpY, tmpZ);

    if (VERBOSE)
        std::cout << "DONE.\n";

    return center;
}

// computes volume of a convex triangulated mesh with centroid "centroid"
double compute_volume(Polyhedron P, Point centroid)
{
    if (VERBOSE)
        std::cout << "> Computing Volume...";

    // Use fixed point (centroid) to generate tetrahedron with
    // each triangle of the hull and sum up volumes.
    double ret = 0.0;
    for (Polyhedron::Facet_iterator iter = P.facets_begin(); iter != P.facets_end(); iter++)
    {
        Facet f = *iter;
        assert(f.is_triangle());
        auto edges = iter->halfedge();
        auto a = edges->vertex()->point();
        edges = edges->next();
        auto b = edges->vertex()->point();
        edges = edges->next();
        auto c = edges->vertex()->point();
        Tetrahedron t = Tetrahedron(centroid, a, b, c);
        ret += abs(t.volume());
    }

    if (VERBOSE)
        std::cout << "DONE.\n";

    return ret;
}

double vector_length(std::vector<double> &v)
{
    double ret = 0.0;
    for (unsigned i = 0; i < v.size(); i++) ret += v[i] * v[i];
    return sqrt(ret);
}

void scale_vector(std::vector<double> &v, double scale)
{
    for (unsigned i = 0; i < v.size(); i++) v[i] = v[i] * scale;
}

// see: "An Iterative Method for Reconstructing Convex Polyhedra from Extended Gaussian Images" - James J. I, Little*
// single step in the optimization algorithm. Equation (3) in the paper explains this.
void optimization_step(
    std::vector<double> A_G,
    std::vector<double> G_L,
    std::vector<double> &L,
    double oldf,
    unsigned k,
    ErrorFunction ef
)
{
    // Calculate next step
    if (VERBOSE)
        std::cout << "> Optimization step...";

    // compute <A, G(L)>
    double GdotA = 0.0;
    for (unsigned i = 0; i < L.size(); i++) GdotA += A_G[i] * G_L[i];

    // compute <A, G(L)> * G(L) - A
    std::vector<double> Stepping;
    Stepping.resize(L.size());
    for (unsigned i = 0; i < L.size(); i++) Stepping[i] = ((GdotA * G_L[i]) - A_G[i]);

    double gamma = 1.0 - GdotA;
    for (unsigned j = 0; j < G_L.size(); j++) L[j] = L[j] + (Stepping[j] * gamma);
}

// see: "An Iterative Method for Reconstructing Convex Polyhedra from Extended Gaussian Images" - James J. I, Little*
int little(std::vector<double> &L, ErrorFunction ef)
{
    unsigned long k = 0;
    double previous = DBL_MAX;
    double current = DBL_MAX;
    Polyhedron P;

    std::cout << "Computing Caustic\n";

    while (true)
    {
        if (VERBOSE)
        {
            std::cout << "> Step: " << k << ":" << '\n';
        }

        std::vector<Kernel::Plane_3> halfspaces = ef.halfspaces(L);
        P = construct_polytope(halfspaces);

        if (VERBOSE)
        {
            std::ofstream out("iteration" + std::to_string(k) + ".off");
            out << P;
        }

        // compute centroid
        Point centroid = compute_centroid(P);

        // translate P(L) to centroid
        if (VERBOSE)
            std::cout << "> Transforming to " << centroid.x() << "," << centroid.y() << "," << centroid.z() << "...";

        Vector translate = Vector(centroid, Point(0, 0, 0));
        for (auto v = P.vertices_begin(); v != P.vertices_end(); v++) v->point() = v->point() + translate;

        if (VERBOSE)
            std::cout << "DONE.\n";

        // Compute V(L)
        double V = compute_volume(P, Point(0, 0, 0));

        // Scale L such that V(L) == 1
        if (VERBOSE)
            std::cout << "> Scaling L...";

        double scale = 1 / pow(V, 0.333);
        scale_vector(L, scale);

        if (VERBOSE)
            std::cout << "DONE.\n";

        // save old value for later comparison
        previous = current;
        current = ef.f(L);

        if (VERBOSE)
        {
            std::cout << "> f(L): " << ef.f(L) << '\n';
            std::cout << "> delta(f): " << (previous - current) << '\n';
        }

        // if we stepped over we revert and wait for reversions to become big enough
        if (previous - current < EPSILON)
        {
            std::cout << "DONE.\n> CONVERGED AT " << current << std::endl;
            break;
        }

        // scale down P to calculate G
        for (auto v = P.vertices_begin(); v != P.vertices_end(); v++)
        {
            Vector v2 = Vector(Point(0, 0, 0), v->point()) * scale;
            v->point() = Point(v2.x(), v2.y(), v2.z());
        }

        // Output the 2D complex to an OFF file.
        if (VERBOSE)
        {
            std::ofstream out(std::string(FILENAME) + "iteration" + std::to_string(k) + "V1.off");
            out << P;
            std::cout << "DONE.\n> Output: iteration" << k << ".off\n";
        }
        else if (k % 100 == 0)
        {
            std::ofstream out(std::string(FILENAME) + "iteration" + std::to_string(k) + "V1.off");
            out << P;
        }

        // Compute G(L)
        std::vector<double> G_L = ef.G(P);
        scale_vector(G_L, 1 / vector_length(G_L));

        std::vector<double> A_G = ef.A_G;

        // if (VERBOSE)
        // {
        //     std::cout << "\nL:   ";
        //     print_vector(L);
        //     std::cout << "\nG(L):";
        //     print_vector(G_L);
        //     std::cout << "\nA_G: ";
        //     print_vector(A_G);
        //     std::cout << '\n';
        // }

        optimization_step(A_G, G_L, L, current, k, ef);
        k++;
    }

    std::cout << "Final Value of f(x): " << ef.f(L) << std::endl;

    // output result
    std::vector<Kernel::Plane_3> halfspaces = ef.halfspaces(L);
    P = construct_polytope(halfspaces);
    std::cout << "> Writing output to final.off\n";
    std::ofstream out(std::string(FILENAME) + "final.off");
    out << P;
    out.close();
    return k;
}

void test_little()
{
    std::vector<double> n1 = { 1, 0, 0 };
    std::vector<double> n2 = { -1, 0, 0 };
    std::vector<double> n3 = { 0, 1, 0 };
    std::vector<double> n4 = { 0, -1, 0 };
    std::vector<double> n5 = { 0, 0, 1 };
    std::vector<double> n6 = { 0, 0, -1 };
    std::vector<std::vector<double>> n = { n1, n2, n3, n4, n5, n6 };
    std::vector<double> A = { 0.5, 0.5, 1, 1, 1, 1 };
    std::vector<double> x = { 1, 1, 1, 1, 1, 1 };
    ErrorFunction ef = ErrorFunction(n, A);
    little(x, ef);

    // vector<double> n1 = {1,0,0};
    // vector<double> n2 = {-1,0,0};
    // vector<double> n3 = {0,1,0};
    // vector<double> n4 = {0,-1,0};
    // vector<double> n5 = {0,0,1};
    // vector<double> n6 = {0,0,-1};
    // vector<double> n7 = {sqrt(0.5),sqrt(0.5),0};
    // vector<vector<double>> n = {n1, n2, n3, n4, n5, n6, n7};
    // vector<double> A = {1, 1+sqrt(0.5), 1, 1+sqrt(0.5), 1, 1, 1};
    // vector<double> x = {1,1,1,1,1,1,1};
    // vector<double> egierror = {0.0, 0.0, 0.0};
    // for (unsigned i = 0; i < n.size(); i++) {
    //     egierror[0] += A[i] * n[i][0];
    //     egierror[1] += A[i] * n[i][1];
    //     egierror[2] += A[i] * n[i][2];
    // }
    // cout << "EGI_Error: (" << egierror[0] << ',' << egierror[1] << ',' << egierror[2] << ")\n";
    // ErrorFunction ef = ErrorFunction(n, A);
    // little(x, ef);
}
