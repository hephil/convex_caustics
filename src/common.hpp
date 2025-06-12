#pragma once

#include <vector>
#include <iostream>
#include <string>

/**
 * @brief Represents a greyscale image with linear brightness values in row-major order in a float pixel array.
 */
struct image_t
{
    /// Width of the input image in pixels.
    uint32_t width;
    /// Height of the input image in pixels.
    uint32_t height;
    /// A flat array of pixel data (e.g. grayscale or luminance values), ordered row-major.
    std::vector<float> pixels;
};

/**
 * @brief Transforms a colour value in sRGB colour space into linear colour space.
 *
 * @param c Colour value in sRGB colour space.
 *
 * @return Colour value in linear colour space.
 */
inline float srgb_to_linear(float c)
{
    const float SRGB_ALPHA = 0.055;
    if (c <= 0.04045)
        return c / 12.92;
    else
        return std::pow((c + SRGB_ALPHA) / (1.0 + SRGB_ALPHA), 2.4);
}

inline void ReleaseAssert(bool c, std::string msg)
{
    if (!c)
    {
        std::cerr << "FATAL: " << msg << std::endl;
        std::abort();
    }
}

#ifdef NDEBUG
#define Assert(c, msg) ;
#else
#define STRINGIFY(s) #s
#define FILELINE(line) __FILE__ ":" STRINGIFY(line) ": "

#ifdef _WIN32
#include <crtdbg.h>
#define Assert(c, msg)                                           \
    {                                                            \
        if (!(c))                                                \
        {                                                        \
            std::cerr << FILELINE(__LINE__) << msg << std::endl; \
            __debugbreak();                                      \
        }                                                        \
    }
#else
#define Assert(c, msg)                                           \
    {                                                            \
        if (!(c))                                                \
        {                                                        \
            std::cerr << FILELINE(__LINE__) << msg << std::endl; \
            std::abort()                                         \
        }                                                        \
    }
#endif
#endif
