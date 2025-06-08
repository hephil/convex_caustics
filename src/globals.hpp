#pragma once

#include <iostream>
#include <string>

#define EPSILON 1e-9
extern bool VERBOSE;
extern unsigned CAUSTIC_SIZE;
extern unsigned PIXEL_SIZE;
extern unsigned DISTANCE;
extern unsigned NUM_GREYVALUES;
extern char *FILENAME;

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
