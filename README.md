# CT Reconstruction - Pure C++ Implementation

Complete implementation of SART, SIRT, and MART algorithms for CT reconstruction in pure C++.

## Features

✅ **Three Reconstruction Algorithms:**
- **SART** (Simultaneous Algebraic Reconstruction Technique) - Best for quality
- **SIRT** (Simultaneous Iterative Reconstruction Technique) - Balanced approach
- **MART** (Multiplicative Algebraic Reconstruction Technique) - Alternative method

✅ **Advanced Features:**
- Fan-beam geometry support
- Negative-log preprocessing
- Multi-threaded processing (OpenMP)
- Command-line interface
- Progress monitoring with time estimates
- Bilinear interpolation
- Non-negativity constraints

✅ **No External Dependencies** (except standard libraries):
- Only requires: libtiff, OpenMP, C++17 compiler

## Quick Start

### 1. Setup and Build

```bash
# Create project directory
cd ~
mkdir -p ct_reconstruction && cd ct_reconstruction

# Copy the source files:
# - ct_reconstruction.cpp (from artifact above)
# - CMakeLists.txt (from artifact above)

# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### 2. Prepare Your Data

```bash
# Create input directory
mkdir -p ../raw/projections

# Copy your .tif projection files to:
# ~/ct_reconstruction/raw/projections/

# Files should be named sequentially (e.g., proj_0001.tif, proj_0002.tif, ...)
```

### 3. Run Reconstruction

```bash
# Basic usage (SART algorithm, 50 iterations)
./ct_recon

# With custom parameters
./ct_recon --algorithm SART --iterations 100 --relaxation 0.7

# SIRT algorithm
./ct_recon --algorithm SIRT --iterations 100 --relaxation 0.5

# MART algorithm
./ct_recon --algorithm MART --iterations 30 --relaxation 0.3

# Reconstruct specific slice range
./ct_recon --z-min 100 --z-max 200

# See all options
./ct_recon --help
```

## Command-Line Options

```
--proj-folder PATH      Input projection folder (default: ./raw/projections)
--output-folder PATH    Output reconstruction folder (default: ./raw/reconstruction_1)
--algorithm ALG         Algorithm: SART, SIRT, or MART (default: SART)
--iterations N          Number of iterations (default: 50)
--relaxation R          Relaxation parameter 0.0-1.0 (default: 0.5)
--z-min N              Minimum slice index (default: 0)
--z-max N              Maximum slice index (default: 999)
--threads N            Number of threads (default: auto)
--no-log               Skip negative-log preprocessing
--help                 Show help message
```

## Algorithm Comparison

| Algorithm | Speed | Quality | Best Use Case |
|-----------|-------|---------|---------------|
| **SART** | Fast | High | General purpose, best quality per iteration |
| **SIRT** | Medium | High | Smoother results, less noise |
| **MART** | Slow | Good | Specific applications, multiplicative updates |

### Recommended Parameters

**SART:**
- Iterations: 50-100
- Relaxation: 0.5-0.7
- Best for: Fast, high-quality reconstruction

**SIRT:**
- Iterations: 100-200
- Relaxation: 0.3-0.5
- Best for: Smooth, low-noise results

**MART:**
- Iterations: 20-50
- Relaxation: 0.1-0.3
- Best for: Specialized applications

## Geometry Configuration

Edit the `GeometryConfig` struct in `ct_reconstruction.cpp`:

```cpp
struct GeometryConfig {
    int num_projections = 360;        // Number of projection angles
    int detector_rows = 1000;         // Detector height in pixels
    int detector_cols = 1000;         // Detector width in pixels
    float sod_mm = 160.0f;           // Source-to-object distance (mm)
    float sdd_mm = 200.0f;           // Source-to-detector distance (mm)
    float pixel_size_mm = 0.048f;    // Detector pixel size (mm)
    float cor_pixel = 518.0f;        // Center of rotation (pixels)
    float voxel_size_mm = 0.0384f;   // Reconstruction voxel size (mm)
    int recon_xy_dim = 1000;         // Reconstruction image size
};
```

## Performance Tips

1. **Use Release Build:**
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```

2. **Optimize Thread Count:**
   ```bash
   # Auto-detect (default)
   ./ct_recon
   
   # Manual setting
   ./ct_recon --threads 8
   ```

3. **Process Subset for Testing:**
   ```bash
   # Test on 10 slices first
   ./ct_recon --z-min 500 --z-max 510 --iterations 20
   ```

4. **Algorithm Selection:**
   - SART: Fastest convergence
   - SIRT: Best quality/smoothness
   - MART: Specialized cases

## Output

Reconstructed slices are saved as:
```
raw/reconstruction_1/ALGORITHM_reconstruction/slice_XXXX.tif
```

Example:
```
raw/reconstruction_1/SART_reconstruction/slice_0000.tif
raw/reconstruction_1/SART_reconstruction/slice_0001.tif
...
```

## Troubleshooting

### Build Errors

**Error: `filesystem` not found**
```bash
# Update compiler or add explicit linking
target_link_libraries(ct_recon stdc++fs)
```

**Error: TIFF library not found**
```bash
# Ubuntu/Debian
sudo apt-get install libtiff-dev

# CentOS/RHEL
sudo yum install libtiff-devel

# macOS
brew install libtiff
```

**Error: OpenMP not found**
```bash
# Usually included with gcc/clang
# For macOS with Homebrew:
brew install libomp
```

### Runtime Issues

**Out of memory:**
- Reduce number of slices with `--z-min` and `--z-max`
- Reduce reconstruction dimensions in config
- Process in batches

**Slow performance:**
- Use `--threads` to adjust thread count
- Reduce iterations for testing
- Use SART instead of SIRT/MART

**Poor quality:**
- Increase iterations
- Adjust relaxation parameter
- Try different algorithms
- Check input data preprocessing

## Example Workflow

```bash
# 1. Quick test (10 slices, 20 iterations)
./ct_recon --algorithm SART --iterations 20 --z-min 500 --z-max 510

# 2. If good, full reconstruction
./ct_recon --algorithm SART --iterations 50

# 3. Try SIRT for smoother results
./ct_recon --algorithm SIRT --iterations 100 --relaxation 0.4

# 4. Compare algorithms
./ct_recon --algorithm MART --iterations 30 --relaxation 0.2
```

## Technical Details

### Forward Projection
- Ray-driven approach
- Bilinear interpolation
- Sub-voxel sampling (0.5 voxel steps)

### Back Projection
- Weighted contribution along rays
- Atomic operations for thread safety
- Normalized by ray weights

### Convergence
- Non-negativity constraint
- Relaxation parameter control
- Projection-by-projection (SART) or simultaneous (SIRT) updates

## License

This is a research/educational implementation. Modify as needed for your application.

## References

- Andersen, A. H., & Kak, A. C. (1984). "Simultaneous algebraic reconstruction technique (SART)"
- Gilbert, P. (1972). "Iterative methods for the three-dimensional reconstruction of an object from projections"
- Gordon, R., Bender, R., & Herman, G. T. (1970). "Algebraic Reconstruction Techniques (ART)"

## Contact

For issues or questions about the implementation, refer to the ASTRA Toolbox documentation for algorithm details.