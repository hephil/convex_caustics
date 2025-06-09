#pragma once

#include "common.hpp"

#include <iostream>
#include <fstream>

/**
 * @brief Reads a PGM (Portable GrayMap) file of type P2 (ASCII encoding).
 *
 * This function reads a grayscale image stored in a plain-text PGM file. Only the P2 format is supported. The image
 * data is stored in an `image_t` structure.
 *
 * @param filename Path to the PGM file to read.
 * @return image_t Struct containing the image width, height, and pixel data.
 */
inline image_t read_pgm_file(const char *filename)
{
    image_t img{
        .width = 0,
        .height = 0,
        .pixels{},
    };

    // Open file
    std::ifstream inputFile(filename);
    ReleaseAssert(inputFile.is_open(), "Unable to open image\n");

    // Read type of PGM
    std::string type;
    inputFile >> type;

    // Plain (ASCII) format
    if (type == "P2")
    {
        float maxValue = 255.0f;
        inputFile >> img.width >> img.height >> maxValue;
        img.pixels.resize(img.width * img.height);
        for (uint32_t i = 0; i < img.width * img.height; i++)
        {
            uint32_t val;
            inputFile >> val;
            img.pixels[i] = srgb_to_linear(val / maxValue);
        }
    }
    else
    {
        std::cout << "Only P2 type pgm files are supported. Provided file has type: " << type << '\n';
    }
    return img;
}
