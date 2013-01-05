cmake_minimum_required(VERSION 2.8)
PROJECT ( ldscpp )
# C++11 and store symbol information
set(CMAKE_CXX_FLAGS "-std=c++11 -g")
set(CMAKE_BUILD_TYPE Release)

find_package(nputils REQUIRED)

include_directories(${NPUTILS_INCLUDE_DIRS})
link_directories(${NPUTILS_LINK_DIRS})


add_subdirectory(unit3d)
add_subdirectory(ang2d)
add_subdirectory(boss)