cmake_minimum_required(VERSION 2.8)
PROJECT ( ldscpp )
set(CMAKE_CXX_FLAGS -std=c++11)
set(CMAKE_BUILD_TYPE Release)

find_package(nputils REQUIRED)

include_directories(${NPUTILS_INCLUDE_DIRS})
link_directories(${NPUTILS_LINK_DIRS})


add_subdirectory(unit3d)