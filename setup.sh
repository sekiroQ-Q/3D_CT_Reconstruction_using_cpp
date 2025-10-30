#!/bin/bash

echo "================================================================"
echo "   CT Reconstruction - Complete Setup Script"
echo "================================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Project directory
PROJECT_DIR="$HOME/ct_reconstruction"

echo "Creating project directory: $PROJECT_DIR"
mkdir -p "$PROJECT_DIR"
cd "$PROJECT_DIR"

# Check for required tools
echo ""
echo "Checking for required tools..."
check_tool() {
    if command -v $1 &> /dev/null; then
        echo -e "${GREEN}✓${NC} $1 found: $(which $1)"
        return 0
    else
        echo -e "${RED}✗${NC} $1 not found"
        return 1
    fi
}

MISSING_TOOLS=0
check_tool cmake || MISSING_TOOLS=1
check_tool g++ || MISSING_TOOLS=1
check_tool make || MISSING_TOOLS=1

# Check for libraries
echo ""
echo "Checking for required libraries..."

# Check for libtiff
if pkg-config --exists libtiff-4 2>/dev/null; then
    echo -e "${GREEN}✓${NC} libtiff found"
else
    echo -e "${YELLOW}⚠${NC} libtiff not found"
    echo "  Installing libtiff..."
    if command -v apt-get &> /dev/null; then
        sudo apt-get update && sudo apt-get install -y libtiff-dev
    elif command -v yum &> /dev/null; then
        sudo yum install -y libtiff-devel
    elif command -v brew &> /dev/null; then
        brew install libtiff
    else
        echo -e "${RED}Error: Unable to install libtiff automatically${NC}"
        MISSING_TOOLS=1
    fi
fi

if [ $MISSING_TOOLS -eq 1 ]; then
    echo ""
    echo -e "${RED}Error: Missing required dependencies${NC}"
    echo "Please install them and run this script again"
    exit 1
fi

echo ""
echo "================================================================"
echo "   Creating source files..."
echo "================================================================"

# Save the full C++ source
cat > ct_reconstruction.cpp << 'CPPEOF'
// Full implementation will be here - copy from the artifact above
// For brevity, indicating where to paste the complete code
#error "Please copy the complete ct_reconstruction.cpp code from the artifact above"
CPPEOF

echo -e "${GREEN}✓${NC} ct_reconstruction.cpp created (needs full code)"

# Save CMakeLists.txt
cat > CMakeLists.txt << 'CMAKEEOF'
cmake_minimum_required(VERSION 3.10)
project(CTReconstruction CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")

find_package(TIFF REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${TIFF_INCLUDE_DIRS})

add_executable(ct_recon ct_reconstruction.cpp)

target_link_libraries(ct_recon
    ${TIFF_LIBRARIES}
    OpenMP::OpenMP_CXX
)

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(ct_recon stdc++fs)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "TIFF: ${TIFF_LIBRARIES}")
message(STATUS "OpenMP: ${OpenMP_FOUND}")
CMAKEEOF

echo -e "${GREEN}✓${NC} CMakeLists.txt created"

# Create directory structure
mkdir -p raw/projections
mkdir -p raw/reconstruction_1

echo -e "${GREEN}✓${NC} Directory structure created"

echo ""
echo "================================================================"
echo "   Building project..."
echo "================================================================"

# Build
mkdir -p build
cd build

echo "Running CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

if [ $? -eq 0 ]; then
    echo ""
    echo "Running Make..."
    make -j$(nproc)
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "================================================================"
        echo -e "${GREEN}   BUILD SUCCESSFUL!${NC}"
        echo "================================================================"
        echo ""
        echo "Executable: $PROJECT_DIR/build/ct_recon"
        echo ""
        echo "Next steps:"
        echo "1. Copy the full source code into ct_reconstruction.cpp"
        echo "2. Rebuild with: cd build && make"
        echo "3. Place .tif files in: $PROJECT_DIR/raw/projections/"
        echo "4. Run: ./ct_recon --help"
        echo ""
        echo "Example commands:"
        echo "  ./ct_recon --algorithm SART --iterations 50"
        echo "  ./ct_recon --algorithm SIRT --iterations 100 --relaxation 0.7"
        echo "  ./ct_recon --algorithm MART --iterations 30 --relaxation 0.3"
        echo ""
    else
        echo -e "${RED}Build failed${NC}"
        exit 1
    fi
else
    echo -e "${RED}CMake configuration failed${NC}"
    exit 1
fi