#include "egi.hpp"

#include <format>
#include <iostream>

void print_egi_stats(discrete_egi_t egi)
{
    std::array<double, 3> egi_error = { 0.0, 0.0, 0.0 };
    egi_error = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < egi.normals.size(); i++)
    {
        egi_error = {
            egi_error[0] - egi.normals[i][0] * egi.areas[i],
            egi_error[1] - egi.normals[i][1] * egi.areas[i],
            egi_error[2] - egi.normals[i][2] * egi.areas[i],
        };
    }
    double accarea = 0.0;
    for (const double w : egi.areas) accarea += w;

    std::cout << "EGI:\n";
    std::cout << std::format("\tError: ({}, {}, {})\n", egi_error[0], egi_error[1], egi_error[2]);
    std::cout << std::format(
        "\tError face: ({}, {}, {})\n",
        egi.normals[egi.normals.size() - 1][0],
        egi.normals[egi.normals.size() - 1][1],
        egi.normals[egi.normals.size() - 1][2]
    );
    std::cout << std::format("\tNumber of facets: {}\n", egi.normals.size() - 2);
    std::cout << std::format(
        "\tArea:\n\t\tfacets: {}\n\t\tbacksize: {}\n\t\tcompensation: {}\n",
        accarea - egi.areas[egi.areas.size() - 1] - egi.areas[egi.areas.size() - 2],
        egi.areas[egi.areas.size() - 2],
        egi.areas[egi.areas.size() - 1]
    );
}

discrete_egi_t make_convex_egi(uint32_t width, uint32_t height, std::span<float> pixels)
{
    const double distance = 1.0; // distance of the projected caustic from the lens
    const double pixel_size = 1.0 / std::max(width, height);
    const double caustic_size = pixel_size * width * pixel_size * height;

    // since we refract all light the total brightness of the projected image is fixed.
    // therefore, we normalise the input image to a brightness value of 1
    double total_brightness = 0;
    size_t num_pixels = std::min(pixels.size(), static_cast<size_t>(width * height));
    size_t num_front_faces = 0u;
    for (size_t i = 0; i < num_pixels; i++)
    {
        total_brightness += pixels[i];
        num_front_faces += pixels[i] > 0.0f ? 1 : 0;
    }

    // the target lense consists of 3 parts:
    // 1. the front face has a distribution of surface normals following the brighntess distribution of the image
    // 2. a compensation surface parallel to the incident light direction to make the lense plausible
    // 3. the backface has an area equal the weighted sum of projected area of the front face normals
    discrete_egi_t egi;
    egi.normals.reserve(num_front_faces + 2u);
    egi.areas.reserve(num_front_faces + 2u);

    double tl_x = -((static_cast<double>(width) - 1) / 2.0) * pixel_size;
    double tl_y = -((static_cast<double>(height) - 1) / 2.0) * pixel_size;
    for (size_t i = 0; i < pixels.size(); i++)
    {
        if (pixels[i] <= 0.0)
            continue;

        uint32_t x = i % width;
        uint32_t y = i / width;

        double nx = tl_x + x * pixel_size;
        double ny = tl_y + y * pixel_size;
        double nz = distance;
        double len = std::sqrt(nx * nx + ny * ny + nz * nz);
        nx /= len;
        ny /= len;
        nz /= len;

        double area = caustic_size * pixels[i] / std::max(1e-20, (std::abs(nz) * total_brightness));

        egi.normals.push_back({ nx, ny, nz });
        egi.areas.push_back(area);
    }

    // flat backside
    egi.normals.push_back({ 0.0, 0.0, -1.0 });
    egi.areas.push_back(caustic_size);

    // Calculate total area of facets and backside, thenadd compensation area based on how lopsided the caustic is
    std::array<double, 3> egi_error = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < egi.normals.size(); i++)
    {
        egi_error = {
            egi_error[0] - egi.normals[i][0] * egi.areas[i],
            egi_error[1] - egi.normals[i][1] * egi.areas[i],
            egi_error[2] - egi.normals[i][2] * egi.areas[i],
        };
    }
    double error_area = std::sqrt(
        egi_error[0] * egi_error[0] + egi_error[1] * egi_error[1] + egi_error[2] * egi_error[2]
    );
    egi.normals.push_back(
        {
            egi_error[0] / error_area,
            egi_error[1] / error_area,
            egi_error[2] / error_area,
        }
    );
    egi.areas.push_back(error_area);

#ifndef NDEBUG
    print_egi_stats(egi);
#endif
    return std::move(egi);
}
