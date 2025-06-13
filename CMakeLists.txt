cmake_minimum_required(VERSION 3.30 FATAL_ERROR)

project (convex_caustics LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftime-trace")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")


set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

find_package(CGAL REQUIRED)

file(GLOB SRC "src/*.inl" "src/*.cpp" "src/*.hpp")
add_executable(convcaust ${SRC})
target_link_libraries(convcaust PRIVATE CGAL::CGAL)

if(LOGGING_ENABLED)
    target_compile_definitions(convcaust PRIVATE LOGGING_ENABLED)
endif()