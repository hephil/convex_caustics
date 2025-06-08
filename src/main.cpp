#include "egi.hpp"
#include "globals.hpp"
#include "little.hpp"

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

bool VERBOSE = false;
unsigned CAUSTIC_SIZE = 100; // 10x10mm
unsigned PIXEL_SIZE = 10;    // 1x1mmm
unsigned DISTANCE = 1000;    // 1000mm
unsigned NUM_GREYVALUES = 64;
char *FILENAME = NULL;

void usage()
{
    cout
        << "usage:\n"
           "convcaust FILENAME.pgm\n output"
           "Options:\n"
           "    -v:    verbose\n"
           "    -d N:  distance between Wall and caustic\n"
           "    -c N:  Area of resulting caustic in mm²\n"
           "    -p N:  Area of a pixel in the projected image in mm²\n"
           "    -g N:  Number of greyvalues the projected image should have\n";
    exit(0);
}

int main(int argc, char **argv)
{
    bool testing = false;

    if (argc < 3)
        usage(); // exits

    FILENAME = argv[argc - 1];

    for (int i = 1; i < argc - 2; i++)
    {
        if (!strcmp(argv[i], "-v"))
        {
            VERBOSE = true;
            continue;
        }
        if (!strcmp(argv[i], "-c"))
        {
            if (i + 2 >= argc)
                usage(); // exits
            long buf = strtol(argv[++i], NULL, 10);
            if (errno || buf <= 0)
                usage(); // exits
            CAUSTIC_SIZE = (unsigned)buf;
            continue;
        }
        if (!strcmp(argv[i], "-p"))
        {
            if (i + 2 >= argc)
                usage(); // exits
            long buf = strtol(argv[++i], NULL, 10);
            if (errno || buf <= 0)
                usage(); // exits
            PIXEL_SIZE = (unsigned)buf;
            continue;
        }
        if (!strcmp(argv[i], "-d"))
        {
            if (i + 2 >= argc)
                usage(); // exits
            long buf = strtol(argv[++i], NULL, 10);
            if (errno || buf <= 0)
                usage(); // exits
            DISTANCE = (unsigned)buf;
            continue;
        }
        if (!strcmp(argv[i], "-g"))
        {
            if (i + 2 >= argc)
                usage(); // exits
            long buf = strtol(argv[++i], NULL, 10);
            if (errno || buf <= 0)
                usage(); // exits
            NUM_GREYVALUES = (unsigned)buf;
            continue;
        }
        if (!strcmp(argv[i], "-t"))
        {
            testing = true;
            continue;
        }
        usage(); // exits
    }

    if (testing)
    {
        test_little();
        return 0;
    }

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

    ImagePGM image = read_pgm(argv[argc - 2]);

    rescale_greyvalues(&image, NUM_GREYVALUES);

    write_pgm(&image);

    // construct EGI (Extended Gaussian Image)
    vector<vector<double>> n;
    vector<double> A;
    construct_egi(&image, n, A, NUM_GREYVALUES);

    // start-values
    vector<double> x;
    for (unsigned i = 0; i < n.size(); i++) x.push_back(1.0);

    // Calculate geometry from EGI
    ErrorFunction ef = ErrorFunction(n, A);
    little(x, ef);

    return 0;
}
