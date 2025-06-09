#include "egi.hpp"
#include "common.hpp"
#include "little.hpp"
#include "pgm.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

char *FILENAME = NULL;

void usage()
{
    cout
        << "usage:\n"
           "convcaust input.pgm output.off\n";
    exit(0);
}

int main(int argc, char **argv)
{

    if (argc < 3)
        usage(); // exits

    FILENAME = argv[argc - 1];

    int i = 0;
    int pos = -1;
    while (argv[argc - 2][i])
    {
        if (argv[argc - 2][i] == '.')
        {
            pos = i;
        }
        i++;
    }
    if (pos == -1)
        usage(); // exits

    string fexts = argv[argc - 2];
    fexts = fexts.substr(pos + 1, string::npos);

    // load the image
    if (fexts != "pgm")
    {
        cout << "Unsupported fileformat: " << fexts << '\n';
        return -1;
    }

    // read input file
    image_t img = read_pgm_file(argv[argc - 2]);

    // construct EGI (Extended Gaussian Image)
    discrete_egi_t egi = make_convex_egi(img.width, img.height, img.pixels);

    // transform format
    vector<vector<double>> n;
    vector<double> A;
    n.resize(egi.normals.size());
    for (size_t i = 0; i < egi.normals.size(); i++)
        n[i] = std::vector<double>({ egi.normals[i][0], egi.normals[i][1], egi.normals[i][2] });

    A.resize(egi.normals.size());
    for (size_t i = 0; i < A.size(); i++) A[i] = egi.normals[i][3];

    // start-values
    vector<double> x;
    for (unsigned i = 0; i < n.size(); i++) x.push_back(1.0);

    // Calculate geometry from EGI
    ErrorFunction ef = ErrorFunction(n, A);
    little(x, ef);

    return 0;
}
