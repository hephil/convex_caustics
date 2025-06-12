#include "egi.hpp"
#include "common.hpp"
#include "little.hpp"
#include "pgm.hpp"

#include <cstdlib>
#include <iostream>

void usage()
{
    std::cout
        << "usage:\n"
           "convcaust input.pgm output.off\n";
    exit(0);
}

int main(int argc, char **argv)
{

    if (argc < 3)
        usage(); // exits

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

    std::string fexts = argv[argc - 2];
    fexts = fexts.substr(pos + 1, std::string::npos);

    // load the image
    if (fexts != "pgm")
    {
        std::cout << "Unsupported fileformat: " << fexts << '\n';
        return -1;
    }

    // read input file
    image_t img = read_pgm_file(argv[argc - 2]);

    // construct EGI (Extended Gaussian Image)
    discrete_egi_t egi = make_convex_egi(img.width, img.height, img.pixels);

    // compute lens
    auto result_mesh = compute_caustic(egi.normals, egi.areas);

    // write lens
    std::ofstream out(argv[argc - 1]);
    if (!out.is_open())
    {
        std::cout << std::format("Could not open output file ({}), storing in out.off instead!", argv[argc - 1]);
        out = std::ofstream("out.off");
    }
    out << result_mesh;

    return 0;
}
