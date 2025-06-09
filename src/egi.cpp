#include "egi.hpp"

#include <iostream>

static const double CAUSTIC_SIZE = 100.0; // 10x10mm
static const double PIXEL_SIZE = 10;      // 1x1mmm
static const double DISTANCE = 1000.0;    // 1000mm

void print_egi_stats(discrete_egi_t egi)
{
    std::array<double, 3> egi_error = { 0.0, 0.0, 0.0 };
    std::cout << "EGI_Error: (" << egi_error[0] << ',' << egi_error[1] << ',' << egi_error[2] << ")\n";
    egi_error = { 0.0, 0.0, 0.0 };
    for (const std::array<double, 4> n : egi.normals)
        egi_error = {
            egi_error[0] - n[0] * n[3],
            egi_error[1] - n[1] * n[3],
            egi_error[1] - n[1] * n[3],
        };
    std::cout << "EGI_Error: (" << egi_error[0] << ',' << egi_error[1] << ',' << egi_error[2] << ")\n";

    std::cout << "EGI:\nNumber of facets: " << egi.normals.size() << '\n';
    double accarea = 0.0;
    for (const std::array<double, 4> n : egi.normals)
    {
        // std::cout << n[0] << ", " << n[1] << ", " << n[2] << ": " << n[3] << '\n';
        accarea += n[3];
    }
    std::cout << "Area of facets: " << accarea - CAUSTIC_SIZE << " backside: " << CAUSTIC_SIZE << '\n';
}

discrete_egi_t make_convex_egi(uint32_t width, uint32_t height, std::span<float> pixels)
{
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

    double tl_x = -((static_cast<double>(width) - 1) / 2.0) * PIXEL_SIZE;
    double tl_y = -((static_cast<double>(height) - 1) / 2.0) * PIXEL_SIZE;
    for (size_t i = 0; i < pixels.size(); i++)
    {
        if (pixels[i] <= 0.0)
            continue;

        uint32_t x = i % width;
        uint32_t y = i / width;

        double nx = tl_x + x * PIXEL_SIZE;
        double ny = tl_y + y * PIXEL_SIZE;
        double nz = DISTANCE;
        double len = std::sqrt(nx * nx + ny * ny + nz * nz);
        nx /= len;
        ny /= len;
        nz /= len;

        double area = CAUSTIC_SIZE * pixels[i] / std::max(1e-20, (std::fabs(nz) * total_brightness));

        egi.normals.push_back({ nx, ny, nz, area });
    }

    // flat backside
    egi.normals.push_back({ 0.0, 0.0, -1.0, CAUSTIC_SIZE });

    // Calculate and add compensation area based on how lopsided the caustic is
    std::array<double, 3> egi_error = { 0.0, 0.0, 0.0 };
    for (const std::array<double, 4> n : egi.normals)
        egi_error = {
            egi_error[0] - n[0] * n[3],
            egi_error[1] - n[1] * n[3],
            egi_error[1] - n[1] * n[3],
        };
    double error_area = std::sqrt(
        egi_error[0] * egi_error[0] + egi_error[1] * egi_error[1] + egi_error[2] * egi_error[2]
    );
    egi.normals.push_back(
        {
            egi_error[0] / error_area,
            egi_error[1] / error_area,
            egi_error[2] / error_area,
            error_area,
        }
    );

    return std::move(egi);
}
