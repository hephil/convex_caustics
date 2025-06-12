#pragma once

#include <vector>
#include <array>
#include <span>
#include <math.h>

/**
 * @brief Represents a discrete convex Extended Gaussian Image (EGI).
 *
 * A discrete EGI is a set of surface normals, each represented as a 4D vector. The first three components (x, y, z)
 * specify the direction of the surface normal. The fourth component is a weight associated with that normal.
 *
 */
struct discrete_egi_t
{
    /// Surface normals; stored as [x, y, z].
    std::vector<std::array<double, 3>> normals;
    /// Surface normal weights
    std::vector<double> areas;
};

/**
 * @brief Converts an image to a convex, discrete EGI (Extended Gaussian Image) representation.
 *
 * This function processes an input image with the given width and height, using its pixel values (e.g. grayscale
 * intensities) to generate surface orientation information stored in a discrete EGI format.
 *
 * @param width Width of the input image in pixels.
 * @param height Height of the input image in pixels.
 * @param pixels A flat array of pixel data (e.g. grayscale or luminance values), ordered row-major.
 *
 * @return A discrete_egi structure containing the resulting convex surface normals.
 *
 * @note The input vector should contain exactly (width * height) elements.
 */
discrete_egi_t make_convex_egi(uint32_t width, uint32_t height, std::span<float> pixels);
