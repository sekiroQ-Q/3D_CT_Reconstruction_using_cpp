cmake_minimum_required(VERSION 3.10)
project(CTReconstruction CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

# Find required packages
find_package(TIFF REQUIRED)
find_package(OpenMP REQUIRED)

# Include directories
include_directories(
    ${TIFF_INCLUDE_DIRS}
)

# Create executable
add_executable(ct_recon
    ct_reconstruction.cpp
)

# Link libraries
target_link_libraries(ct_recon
    ${TIFF_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# For older compilers that need explicit filesystem lib
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(ct_recon stdc++fs)
endif()

# Installation
install(TARGETS ct_recon DESTINATION bin)

# Print build information
message(STATUS "========================================")
message(STATUS "CT Reconstruction Build Configuration")
message(STATUS "========================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "TIFF Found: ${TIFF_FOUND}")
message(STATUS "TIFF Libraries: ${TIFF_LIBRARIES}")
message(STATUS "OpenMP Found: ${OpenMP_FOUND}")
if(OpenMP_FOUND)
    message(STATUS "OpenMP Version: ${OpenMP_CXX_VERSION}")
endif()
message(STATUS "========================================")

# Create directory structure
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/../raw/projections)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/../raw/reconstruction_1)