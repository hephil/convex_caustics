#pragma once

#include <vector>
#include <math.h>

class ImagePGM
{

public:
    std::vector<std::vector<int>> data;
    int width;
    int height;
    int maxValue;

    ImagePGM(std::vector<std::vector<int>>, int, int, int);
};

ImagePGM read_pgm(char *);
void write_pgm(ImagePGM *);
void rescale_greyvalues(ImagePGM *, int);
void construct_egi(ImagePGM *, std::vector<std::vector<double>> &, std::vector<double> &, unsigned numValues);
void normalize3(std::vector<double> &);
