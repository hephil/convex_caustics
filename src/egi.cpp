#include "globals.hpp"
#include "egi.hpp"
#include <fstream>

ImagePGM::ImagePGM(std::vector<std::vector<int>> image, int width, int height, int maxValue)
{
    this->data = image;
    this->width = width;
    this->height = height;
    this->maxValue = maxValue;
}

// Reads a .pgm file of type P2.
ImagePGM read_pgm(char *filename)
{
    std::string type;
    int height;
    int width;
    int maxValue;
    std::vector<std::vector<int>> image;

    // open file
    std::ifstream inputFile(filename);
    ReleaseAssert(inputFile.is_open(), "Unable to open image\n");

    // read type of pgm
    inputFile >> type;

    // plain format
    if (type == "P2")
    {
        inputFile >> width >> height >> maxValue;

        image.resize(height);

        for (int i = 0; i < width; i++)
        {
            image[i].resize(width);

            for (int j = 0; j < height; j++)
            {
                inputFile >> image[i][j];
            }
        }

        return ImagePGM(image, width, height, maxValue);
    }
    // invalid type
    else
    {
        std::cout << "Invalid .pgm type: " << type << '\n';
        exit(0);
        // exits
    }
}

// just to test the reading and rescaling, thus output hardcoded. Copied from Didyk.
void write_pgm(ImagePGM *image)
{

    std::ofstream outputFile("out.pgm");

    if (!outputFile.is_open())
    {
        std::cout << "Unable to open OutputImageFile" << std::endl;
        return;
    }

    outputFile << "P2" << std::endl << image->width << " " << image->height << std::endl << image->maxValue << std::endl;

    for (int i = 0; i < image->height; i++)
    {
        for (int j = 0; j < image->width; j++)
        {
            outputFile << image->data[i][j] << " ";
        }
        outputFile << std::endl;
    }
}

// Samples the greyvalues of the read image to numValues different greyscales
void rescale_greyvalues(ImagePGM *image, int numValues)
{
    float size = ((float)image->maxValue + 1) / (float)numValues;
    for (int x = 0; x < image->width; x++)
    {
        for (int y = 0; y < image->height; y++)
        {
            image->data[x][y] = image->data[x][y] / size;
        }
    }
    image->maxValue = numValues - 1;
}

// Computes the surface normals of the facets. std::vectors are in format (x,y,z,b),
// where x is along the width of the image, y along the height, z is the distance to the caustic, and b is the fraction
// of the overall brightness of the pixel. the center of the coordinate system is the caustic.
void construct_egi(ImagePGM *image, std::vector<std::vector<double>> &normals, std::vector<double> &area, unsigned numValues)
{

    int width = image->width;
    int height = image->height;

    /*
    normals.resize(width * height + 1);
    normals[width * height].resize(3);
    area.resize(width * height + 1);
    */

    double brightness = 0;

    // coordinates of top-left pixel
    double tl_x = -((double)(width - 1) / 2.0) * PIXEL_SIZE;
    double tl_y = -((double)(height - 1) / 2.0) * PIXEL_SIZE;

    // compute overall brightness and resize std::vectors
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            // normals[y * width + x].resize(3);
            double tmp = (double)image->data[x][y] / (double)numValues;
            brightness += pow(tmp, 2);
        }
    }

    int count = 0;
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if (!image->data[x][y])
                continue;
            double x0 = tl_x + x * PIXEL_SIZE;
            double y0 = tl_y + y * PIXEL_SIZE;
            double z0 = DISTANCE;
            std::vector<double> normal = { x0, y0, z0 };
            normalize3(normal);
            normals.push_back(normal);

            // compute area of facet
            // A = (p * CA_SI) / (normal * (0,0,1))
            double tmpp = (double)image->data[x][y] / (double)numValues;
            double p = pow(tmpp, 2) / brightness;
            double a = (CAUSTIC_SIZE * p) / fabs(normals[count][2]);
            area.push_back(a);
            count++;
        }
    }
    // backside
    double bx = 0.0;
    double by = 0.0;
    double bz = -1.0;
    std::vector<double> backside_normal = { bx, by, bz };
    normals.push_back(backside_normal);
    area.push_back(CAUSTIC_SIZE);

    // calculate egi error
    std::vector<double> egi_error = { 0.0, 0.0, 0.0 };
    for (unsigned i = 0; i < normals.size(); i++)
    {
        egi_error[0] -= area[i] * normals[i][0];
        egi_error[1] -= area[i] * normals[i][1];
        egi_error[2] -= area[i] * normals[i][2];
    }

    if (VERBOSE)
        std::cout << "EGI_Error: (" << egi_error[0] << ',' << egi_error[1] << ',' << egi_error[2] << ")\n";

    // fix egi error to make egi feasible
    double error_area = sqrt(egi_error[0] * egi_error[0] + egi_error[1] * egi_error[1] + egi_error[2] * egi_error[2]);
    area.push_back(error_area);
    normalize3(egi_error);
    normals.push_back(egi_error);

    if (VERBOSE)
    {
        egi_error = { 0.0, 0.0, 0.0 };
        for (unsigned i = 0; i < normals.size(); i++)
        {
            egi_error[0] += area[i] * normals[i][0];
            egi_error[1] += area[i] * normals[i][1];
            egi_error[2] += area[i] * normals[i][2];
        }
        std::cout << "EGI_Error: (" << egi_error[0] << ',' << egi_error[1] << ',' << egi_error[2] << ")\n";

        std::ofstream out("picture.egi");
        out << "EGI:\nNumber of facets: " << count << '\n';
        double accarea = 0.0;
        for (int i = 0; i < normals.size(); i++)
        {
            out << normals[i][0] << ", " << normals[i][1] << ", " << normals[i][2] << ": " << area[i] << '\n';
            accarea += area[i];
        }
        std::cout << "Area of facets: " << accarea - CAUSTIC_SIZE << " backside: " << CAUSTIC_SIZE << '\n';

        out.close();
    }
}

void normalize3(std::vector<double> &vec3)
{
    double x = vec3[0];
    double y = vec3[1];
    double z = vec3[2];
    double len = sqrt(x * x + y * y + z * z);

    if (len != 0)
    {
        vec3[0] = x / len;
        vec3[1] = y / len;
        vec3[2] = z / len;
    }
}
