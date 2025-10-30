#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <tiffio.h>
#include <omp.h>
#include <iomanip>
#include <chrono>

namespace fs = std::filesystem;

// Configuration structure
struct GeometryConfig {
    int num_projections = 360;
    int detector_rows = 1000;
    int detector_cols = 1000;
    float sod_mm = 160.0f;           // Source to object distance
    float sdd_mm = 200.0f;           // Source to detector distance
    float pixel_size_mm = 0.048f;    // Detector pixel size
    float cor_pixel = 518.0f;        // Center of rotation in pixels
    float v_center_pixel = 500.0f;   // Vertical center
    int recon_z_min = 0;
    int recon_z_max = 999;
    int recon_xy_dim = 1000;
    int num_iterations = 50;
    float voxel_size_mm = 0.0384f;
    std::string algorithm = "SART";   // SART, SIRT, or MART
    float relaxation = 0.5f;          // Relaxation parameter (0.0-1.0)
    bool apply_log = true;            // Apply negative log preprocessing
    int num_threads = 0;              // 0 = auto-detect
};

// Timer utility
class Timer {
    std::chrono::high_resolution_clock::time_point start_time;
public:
    Timer() : start_time(std::chrono::high_resolution_clock::now()) {}
    
    double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start_time).count();
    }
    
    void reset() {
        start_time = std::chrono::high_resolution_clock::now();
    }
};

// TIFF reading function
bool readTIFF(const std::string& filename, std::vector<float>& data, int& width, int& height) {
    TIFF* tif = TIFFOpen(filename.c_str(), "r");
    if (!tif) {
        std::cerr << "Error: Cannot open TIFF file: " << filename << std::endl;
        return false;
    }
    
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    
    data.resize(width * height);
    
    uint32* raster = (uint32*)_TIFFmalloc(width * height * sizeof(uint32));
    if (raster != nullptr) {
        if (TIFFReadRGBAImage(tif, width, height, raster, 0)) {
            for (int i = 0; i < width * height; i++) {
                uint32 pixel = raster[i];
                float gray = (TIFFGetR(pixel) + TIFFGetG(pixel) + TIFFGetB(pixel)) / 3.0f;
                data[i] = gray;
            }
        }
        _TIFFfree(raster);
    }
    
    TIFFClose(tif);
    return true;
}

// TIFF writing function
bool writeTIFF(const std::string& filename, const float* data, int width, int height) {
    TIFF* tif = TIFFOpen(filename.c_str(), "w");
    if (!tif) {
        std::cerr << "Error: Cannot create TIFF file: " << filename << std::endl;
        return false;
    }
    
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    
    for (int row = 0; row < height; row++) {
        TIFFWriteScanline(tif, (void*)&data[row * width], row, 0);
    }
    
    TIFFClose(tif);
    return true;
}

// Load all projections
bool loadProjections(const std::string& path, 
                     const GeometryConfig& config,
                     std::vector<float>& projections_data) {
    
    std::cout << "Loading projections from: " << path << std::endl;
    Timer timer;
    
    std::vector<std::string> file_list;
    for (const auto& entry : fs::directory_iterator(path)) {
        auto ext = entry.path().extension().string();
        if (ext == ".tif" || ext == ".tiff") {
            file_list.push_back(entry.path().string());
        }
    }
    
    std::sort(file_list.begin(), file_list.end());
    
    if (file_list.empty()) {
        std::cerr << "Error: No .tif files found in '" << path << "'." << std::endl;
        return false;
    }
    
    int num_projections = file_list.size();
    std::cout << "Found " << num_projections << " projection files." << std::endl;
    
    // Allocate memory: (detector_rows, num_projections, detector_cols)
    projections_data.resize(config.detector_rows * num_projections * config.detector_cols);
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_projections; i++) {
        if (i % 50 == 0) {
            #pragma omp critical
            std::cout << "  Loading file " << (i + 1) << " of " << num_projections << "..." << std::endl;
        }
        
        std::vector<float> img_data;
        int width, height;
        
        if (!readTIFF(file_list[i], img_data, width, height)) {
            continue;
        }
        
        if (width != config.detector_cols || height != config.detector_rows) {
            #pragma omp critical
            std::cerr << "Error: Image " << file_list[i] << " has wrong dimensions." << std::endl;
            continue;
        }
        
        // Copy to projections array: projections[row][proj][col]
        for (int row = 0; row < config.detector_rows; row++) {
            for (int col = 0; col < config.detector_cols; col++) {
                int idx = row * num_projections * config.detector_cols + i * config.detector_cols + col;
                projections_data[idx] = img_data[row * width + col];
            }
        }
    }
    
    std::cout << "Successfully loaded " << num_projections << " projections in " 
              << timer.elapsed() << " seconds." << std::endl;
    
    // Preprocessing: Negative-log transform
    if (config.apply_log) {
        std::cout << "\nApplying negative-log preprocessing..." << std::endl;
        
        float i_max = *std::max_element(projections_data.begin(), projections_data.end());
        float i_min = *std::min_element(projections_data.begin(), projections_data.end());
        
        std::cout << "  Raw data range: [" << i_min << ", " << i_max << "]" << std::endl;
        
        if (i_max < 10.0f) {
            std::cout << "  WARNING: Data appears already preprocessed. Skipping log transform." << std::endl;
        } else {
            // Standard log transform
            #pragma omp parallel for
            for (size_t i = 0; i < projections_data.size(); i++) {
                float val = projections_data[i] / i_max;
                val = std::max(1e-6f, std::min(1.0f, val));
                projections_data[i] = -std::log(val);
            }
            
            float new_min = *std::min_element(projections_data.begin(), projections_data.end());
            float new_max = *std::max_element(projections_data.begin(), projections_data.end());
            std::cout << "  After log transform: [" << new_min << ", " << new_max << "]" << std::endl;
        }
    }
    
    return true;
}

// Forward projection for a single ray (fan-beam geometry)
float forwardProjectRay(const std::vector<float>& volume,
                       int vol_size,
                       float source_x, float source_y,
                       float det_x, float det_y) {
    
    float dx = det_x - source_x;
    float dy = det_y - source_y;
    float ray_length = std::sqrt(dx * dx + dy * dy);
    float step_size = 0.5f;  // Sub-voxel sampling
    int num_steps = static_cast<int>(ray_length / step_size) + 1;
    
    float sum = 0.0f;
    float weight_sum = 0.0f;
    
    for (int step = 0; step < num_steps; step++) {
        float t = static_cast<float>(step) / (num_steps - 1);
        float x = source_x + t * dx;
        float y = source_y + t * dy;
        
        // Convert to voxel coordinates (centered)
        float vx = x + vol_size / 2.0f;
        float vy = y + vol_size / 2.0f;
        
        // Bilinear interpolation
        int x0 = static_cast<int>(std::floor(vx));
        int y0 = static_cast<int>(std::floor(vy));
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        
        if (x0 >= 0 && x1 < vol_size && y0 >= 0 && y1 < vol_size) {
            float fx = vx - x0;
            float fy = vy - y0;
            
            float val = (1 - fx) * (1 - fy) * volume[y0 * vol_size + x0] +
                       fx * (1 - fy) * volume[y0 * vol_size + x1] +
                       (1 - fx) * fy * volume[y1 * vol_size + x0] +
                       fx * fy * volume[y1 * vol_size + x1];
            
            sum += val * step_size;
            weight_sum += step_size;
        }
    }
    
    return (weight_sum > 0) ? sum : 0.0f;
}

// Back-projection for a single ray (fan-beam geometry)
void backprojectRay(std::vector<float>& volume,
                   std::vector<float>& weights,
                   int vol_size,
                   float value,
                   float source_x, float source_y,
                   float det_x, float det_y) {
    
    float dx = det_x - source_x;
    float dy = det_y - source_y;
    float ray_length = std::sqrt(dx * dx + dy * dy);
    float step_size = 0.5f;
    int num_steps = static_cast<int>(ray_length / step_size) + 1;
    
    for (int step = 0; step < num_steps; step++) {
        float t = static_cast<float>(step) / (num_steps - 1);
        float x = source_x + t * dx;
        float y = source_y + t * dy;
        
        float vx = x + vol_size / 2.0f;
        float vy = y + vol_size / 2.0f;
        
        int x0 = static_cast<int>(std::floor(vx));
        int y0 = static_cast<int>(std::floor(vy));
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        
        if (x0 >= 0 && x1 < vol_size && y0 >= 0 && y1 < vol_size) {
            float fx = vx - x0;
            float fy = vy - y0;
            
            float w00 = (1 - fx) * (1 - fy) * step_size;
            float w10 = fx * (1 - fy) * step_size;
            float w01 = (1 - fx) * fy * step_size;
            float w11 = fx * fy * step_size;
            
            #pragma omp atomic
            volume[y0 * vol_size + x0] += value * w00;
            #pragma omp atomic
            volume[y0 * vol_size + x1] += value * w10;
            #pragma omp atomic
            volume[y1 * vol_size + x0] += value * w01;
            #pragma omp atomic
            volume[y1 * vol_size + x1] += value * w11;
            
            #pragma omp atomic
            weights[y0 * vol_size + x0] += w00;
            #pragma omp atomic
            weights[y0 * vol_size + x1] += w10;
            #pragma omp atomic
            weights[y1 * vol_size + x0] += w01;
            #pragma omp atomic
            weights[y1 * vol_size + x1] += w11;
        }
    }
}

// SART reconstruction for single slice
void reconstructSliceSART(const std::vector<float>& sinogram,
                         std::vector<float>& reconstruction,
                         const GeometryConfig& config,
                         int slice_idx) {
    
    int vol_size = config.recon_xy_dim;
    int num_det = config.detector_cols;
    int num_proj = config.num_projections;
    
    float sod = config.sod_mm / config.voxel_size_mm;
    float sdd = config.sdd_mm / config.voxel_size_mm;
    float det_spacing = config.pixel_size_mm / config.voxel_size_mm;
    float det_offset = (config.cor_pixel - num_det / 2.0f) * det_spacing;
    
    reconstruction.assign(vol_size * vol_size, 0.0f);
    
    // SART iterations
    for (int iter = 0; iter < config.num_iterations; iter++) {
        Timer iter_timer;
        
        // Process each projection sequentially (SART updates after each projection)
        for (int proj = 0; proj < num_proj; proj++) {
            std::vector<float> correction(vol_size * vol_size, 0.0f);
            std::vector<float> weights(vol_size * vol_size, 1e-6f);
            
            float angle = (2.0f * M_PI * proj) / num_proj;
            float cos_a = std::cos(angle);
            float sin_a = std::sin(angle);
            
            // Source position
            float source_x = -sod * cos_a;
            float source_y = -sod * sin_a;
            
            // Detector center
            float det_center_x = (sdd - sod) * cos_a;
            float det_center_y = (sdd - sod) * sin_a;
            
            // Detector direction (perpendicular to source-detector)
            float det_dir_x = -sin_a;
            float det_dir_y = cos_a;
            
            #pragma omp parallel for
            for (int det = 0; det < num_det; det++) {
                float det_pos = (det - num_det / 2.0f) * det_spacing + det_offset;
                float det_x = det_center_x + det_pos * det_dir_x;
                float det_y = det_center_y + det_pos * det_dir_y;
                
                // Forward project
                float fp = forwardProjectRay(reconstruction, vol_size, 
                                            source_x, source_y, det_x, det_y);
                
                // Compute correction
                float measured = sinogram[proj * num_det + det];
                float diff = measured - fp;
                
                // Back-project correction
                backprojectRay(correction, weights, vol_size, diff,
                              source_x, source_y, det_x, det_y);
            }
            
            // Apply correction for this projection
            #pragma omp parallel for
            for (int i = 0; i < vol_size * vol_size; i++) {
                if (weights[i] > 1e-6f) {
                    reconstruction[i] += config.relaxation * correction[i] / weights[i];
                    reconstruction[i] = std::max(0.0f, reconstruction[i]);  // Non-negativity
                }
            }
        }
        
        if (iter % 10 == 0 || iter == config.num_iterations - 1) {
            float mean = std::accumulate(reconstruction.begin(), reconstruction.end(), 0.0f) / reconstruction.size();
            float min_val = *std::min_element(reconstruction.begin(), reconstruction.end());
            float max_val = *std::max_element(reconstruction.begin(), reconstruction.end());
            
            std::cout << "    Iteration " << std::setw(3) << (iter + 1) << "/" << config.num_iterations 
                     << " | Time: " << std::setw(6) << std::fixed << std::setprecision(2) << iter_timer.elapsed() << "s"
                     << " | Mean: " << std::setw(8) << std::setprecision(4) << mean
                     << " | Range: [" << std::setprecision(4) << min_val << ", " << max_val << "]" << std::endl;
        }
    }
}

// SIRT reconstruction for single slice (simultaneous update)
void reconstructSliceSIRT(const std::vector<float>& sinogram,
                         std::vector<float>& reconstruction,
                         const GeometryConfig& config,
                         int slice_idx) {
    
    int vol_size = config.recon_xy_dim;
    int num_det = config.detector_cols;
    int num_proj = config.num_projections;
    
    float sod = config.sod_mm / config.voxel_size_mm;
    float sdd = config.sdd_mm / config.voxel_size_mm;
    float det_spacing = config.pixel_size_mm / config.voxel_size_mm;
    float det_offset = (config.cor_pixel - num_det / 2.0f) * det_spacing;
    
    reconstruction.assign(vol_size * vol_size, 0.0f);
    
    // SIRT iterations
    for (int iter = 0; iter < config.num_iterations; iter++) {
        Timer iter_timer;
        
        std::vector<float> correction(vol_size * vol_size, 0.0f);
        std::vector<float> weights(vol_size * vol_size, 1e-6f);
        
        // Process all projections (SIRT updates after all projections)
        #pragma omp parallel for
        for (int proj = 0; proj < num_proj; proj++) {
            float angle = (2.0f * M_PI * proj) / num_proj;
            float cos_a = std::cos(angle);
            float sin_a = std::sin(angle);
            
            float source_x = -sod * cos_a;
            float source_y = -sod * sin_a;
            
            float det_center_x = (sdd - sod) * cos_a;
            float det_center_y = (sdd - sod) * sin_a;
            
            float det_dir_x = -sin_a;
            float det_dir_y = cos_a;
            
            for (int det = 0; det < num_det; det++) {
                float det_pos = (det - num_det / 2.0f) * det_spacing + det_offset;
                float det_x = det_center_x + det_pos * det_dir_x;
                float det_y = det_center_y + det_pos * det_dir_y;
                
                float fp = forwardProjectRay(reconstruction, vol_size, 
                                            source_x, source_y, det_x, det_y);
                
                float measured = sinogram[proj * num_det + det];
                float diff = measured - fp;
                
                backprojectRay(correction, weights, vol_size, diff,
                              source_x, source_y, det_x, det_y);
            }
        }
        
        // Apply correction after all projections
        #pragma omp parallel for
        for (int i = 0; i < vol_size * vol_size; i++) {
            if (weights[i] > 1e-6f) {
                reconstruction[i] += config.relaxation * correction[i] / weights[i];
                reconstruction[i] = std::max(0.0f, reconstruction[i]);
            }
        }
        
        if (iter % 10 == 0 || iter == config.num_iterations - 1) {
            float mean = std::accumulate(reconstruction.begin(), reconstruction.end(), 0.0f) / reconstruction.size();
            float min_val = *std::min_element(reconstruction.begin(), reconstruction.end());
            float max_val = *std::max_element(reconstruction.begin(), reconstruction.end());
            
            std::cout << "    Iteration " << std::setw(3) << (iter + 1) << "/" << config.num_iterations 
                     << " | Time: " << std::setw(6) << std::fixed << std::setprecision(2) << iter_timer.elapsed() << "s"
                     << " | Mean: " << std::setw(8) << std::setprecision(4) << mean
                     << " | Range: [" << std::setprecision(4) << min_val << ", " << max_val << "]" << std::endl;
        }
    }
}

// MART reconstruction for single slice (multiplicative update)
void reconstructSliceMART(const std::vector<float>& sinogram,
                         std::vector<float>& reconstruction,
                         const GeometryConfig& config,
                         int slice_idx) {
    
    int vol_size = config.recon_xy_dim;
    int num_det = config.detector_cols;
    int num_proj = config.num_projections;
    
    float sod = config.sod_mm / config.voxel_size_mm;
    float sdd = config.sdd_mm / config.voxel_size_mm;
    float det_spacing = config.pixel_size_mm / config.voxel_size_mm;
    float det_offset = (config.cor_pixel - num_det / 2.0f) * det_spacing;
    
    // Initialize with small positive value
    reconstruction.assign(vol_size * vol_size, 0.1f);
    
    // MART iterations
    for (int iter = 0; iter < config.num_iterations; iter++) {
        Timer iter_timer;
        
        // Process each projection
        for (int proj = 0; proj < num_proj; proj++) {
            float angle = (2.0f * M_PI * proj) / num_proj;
            float cos_a = std::cos(angle);
            float sin_a = std::sin(angle);
            
            float source_x = -sod * cos_a;
            float source_y = -sod * sin_a;
            
            float det_center_x = (sdd - sod) * cos_a;
            float det_center_y = (sdd - sod) * sin_a;
            
            float det_dir_x = -sin_a;
            float det_dir_y = cos_a;
            
            #pragma omp parallel for
            for (int det = 0; det < num_det; det++) {
                float det_pos = (det - num_det / 2.0f) * det_spacing + det_offset;
                float det_x = det_center_x + det_pos * det_dir_x;
                float det_y = det_center_y + det_pos * det_dir_y;
                
                float fp = forwardProjectRay(reconstruction, vol_size, 
                                            source_x, source_y, det_x, det_y);
                
                float measured = sinogram[proj * num_det + det];
                
                // Multiplicative correction factor
                float correction_factor = 1.0f;
                if (fp > 1e-6f && measured > 1e-6f) {
                    correction_factor = std::pow(measured / fp, config.relaxation);
                }
                
                // Apply multiplicative correction along the ray
                float dx = det_x - source_x;
                float dy = det_y - source_y;
                float ray_length = std::sqrt(dx * dx + dy * dy);
                float step_size = 0.5f;
                int num_steps = static_cast<int>(ray_length / step_size) + 1;
                
                for (int step = 0; step < num_steps; step++) {
                    float t = static_cast<float>(step) / (num_steps - 1);
                    float x = source_x + t * dx;
                    float y = source_y + t * dy;
                    
                    float vx = x + vol_size / 2.0f;
                    float vy = y + vol_size / 2.0f;
                    
                    int x0 = static_cast<int>(std::floor(vx));
                    int y0 = static_cast<int>(std::floor(vy));
                    int x1 = x0 + 1;
                    int y1 = y0 + 1;
                    
                    if (x0 >= 0 && x1 < vol_size && y0 >= 0 && y1 < vol_size) {
                        float fx = vx - x0;
                        float fy = vy - y0;
                        
                        float w00 = (1 - fx) * (1 - fy);
                        float w10 = fx * (1 - fy);
                        float w01 = (1 - fx) * fy;
                        float w11 = fx * fy;
                        
                        // Multiplicative update with atomic operations
                        #pragma omp critical
                        {
                            int idx00 = y0 * vol_size + x0;
                            int idx10 = y0 * vol_size + x1;
                            int idx01 = y1 * vol_size + x0;
                            int idx11 = y1 * vol_size + x1;
                            
                            reconstruction[idx00] *= std::pow(correction_factor, w00 * step_size);
                            reconstruction[idx10] *= std::pow(correction_factor, w10 * step_size);
                            reconstruction[idx01] *= std::pow(correction_factor, w01 * step_size);
                            reconstruction[idx11] *= std::pow(correction_factor, w11 * step_size);
                        }
                    }
                }
            }
        }
        
        if (iter % 10 == 0 || iter == config.num_iterations - 1) {
            float mean = std::accumulate(reconstruction.begin(), reconstruction.end(), 0.0f) / reconstruction.size();
            float min_val = *std::min_element(reconstruction.begin(), reconstruction.end());
            float max_val = *std::max_element(reconstruction.begin(), reconstruction.end());
            
            std::cout << "    Iteration " << std::setw(3) << (iter + 1) << "/" << config.num_iterations 
                     << " | Time: " << std::setw(6) << std::fixed << std::setprecision(2) << iter_timer.elapsed() << "s"
                     << " | Mean: " << std::setw(8) << std::setprecision(4) << mean
                     << " | Range: [" << std::setprecision(4) << min_val << ", " << max_val << "]" << std::endl;
        }
    }
}

// Main reconstruction function
bool reconstructSliceBySlice(const std::vector<float>& projections_data,
                            const GeometryConfig& config,
                            const std::string& output_folder) {
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Starting Slice-by-Slice CT Reconstruction" << std::endl;
    std::cout << "Algorithm: " << config.algorithm << std::endl;
    std::cout << "Iterations: " << config.num_iterations << std::endl;
    std::cout << "Relaxation: " << config.relaxation << std::endl;
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
    std::cout << std::string(70, '=') << "\n" << std::endl;
    
    int num_slices = config.recon_z_max - config.recon_z_min + 1;
    fs::create_directories(output_folder);
    
    Timer total_timer;
    
    for (int slice_idx = config.recon_z_min; slice_idx <= config.recon_z_max; slice_idx++) {
        std::cout << "\nSlice " << slice_idx << " (" << (slice_idx - config.recon_z_min + 1) 
                  << "/" << num_slices << ")" << std::endl;
        
        Timer slice_timer;
        
        // Extract sinogram for this slice
        std::vector<float> sinogram(config.num_projections * config.detector_cols);
        for (int proj = 0; proj < config.num_projections; proj++) {
            for (int col = 0; col < config.detector_cols; col++) {
                int idx = slice_idx * config.num_projections * config.detector_cols + 
                         proj * config.detector_cols + col;
                sinogram[proj * config.detector_cols + col] = projections_data[idx];
            }
        }
        
        // Reconstruct this slice using selected algorithm
        std::vector<float> reconstruction;
        
        if (config.algorithm == "SART") {
            reconstructSliceSART(sinogram, reconstruction, config, slice_idx);
        } else if (config.algorithm == "SIRT") {
            reconstructSliceSIRT(sinogram, reconstruction, config, slice_idx);
        } else if (config.algorithm == "MART") {
            reconstructSliceMART(sinogram, reconstruction, config, slice_idx);
        } else {
            std::cerr << "Unknown algorithm: " << config.algorithm << std::endl;
            return false;
        }
        
        // Save the slice
        char filename[256];
        snprintf(filename, sizeof(filename), "slice_%04d.tif", slice_idx);
        std::string save_path = output_folder + "/" + filename;
        
        if (!writeTIFF(save_path, reconstruction.data(), config.recon_xy_dim, config.recon_xy_dim)) {
            std::cerr << "Error: Failed to save slice " << slice_idx << std::endl;
            return false;
        }
        
        float min_val = *std::min_element(reconstruction.begin(), reconstruction.end());
        float max_val = *std::max_element(reconstruction.begin(), reconstruction.end());
        float mean_val = std::accumulate(reconstruction.begin(), reconstruction.end(), 0.0f) / reconstruction.size();
        
        std::cout << "  Saved: " << save_path << std::endl;
        std::cout << "  Statistics: min=" << min_val << ", max=" << max_val 
                  << ", mean=" << mean_val << std::endl;
        std::cout << "  Slice time: " << slice_timer.elapsed() << " seconds" << std::endl;
        
        // Estimate remaining time
        int completed = slice_idx - config.recon_z_min + 1;
        float avg_time = total_timer.elapsed() / completed;
        int remaining = num_slices - completed;
        float est_remaining = avg_time * remaining;
        
        std::cout << "  Estimated remaining time: " << std::fixed << std::setprecision(1) 
                  << est_remaining << " seconds (" << (est_remaining / 60.0) << " minutes)" << std::endl;
    }
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Reconstruction complete!" << std::endl;
    std::cout << "Total time: " << std::fixed << std::setprecision(2) 
              << total_timer.elapsed() << " seconds (" 
              << (total_timer.elapsed() / 60.0) << " minutes)" << std::endl;
    std::cout << num_slices << " slices saved to: " << output_folder << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    return true;
}

// Print usage information
void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]" << std::endl;
    std::cout << "\nOptions:" << std::endl;
    std::cout << "  --proj-folder PATH      Input projection folder (default: ./raw/projections)" << std::endl;
    std::cout << "  --output-folder PATH    Output reconstruction folder (default: ./raw/reconstruction_1)" << std::endl;
    std::cout << "  --algorithm ALG         Algorithm: SART, SIRT, or MART (default: SART)" << std::endl;
    std::cout << "  --iterations N          Number of iterations (default: 50)" << std::endl;
    std::cout << "  --relaxation R          Relaxation parameter 0.0-1.0 (default: 0.5)" << std::endl;
    std::cout << "  --z-min N               Minimum slice index (default: 0)" << std::endl;
    std::cout << "  --z-max N               Maximum slice index (default: 999)" << std::endl;
    std::cout << "  --threads N             Number of threads (default: auto)" << std::endl;
    std::cout << "  --no-log                Skip negative-log preprocessing" << std::endl;
    std::cout << "  --help                  Show this help message" << std::endl;
    std::cout << "\nExample:" << std::endl;
    std::cout << "  " << program_name << " --algorithm SART --iterations 100 --relaxation 0.7" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "            CT Reconstruction Tool - Pure C++ Implementation" << std::endl;
    std::cout << "                   SART / SIRT / MART Algorithms" << std::endl;
    std::cout << std::string(70, '=') << "\n" << std::endl;
    
    // Default configuration
    GeometryConfig config;
    std::string proj_folder = "./raw/projections";
    std::string recon_folder = "./raw/reconstruction_1";
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--proj-folder" && i + 1 < argc) {
            proj_folder = argv[++i];
        } else if (arg == "--output-folder" && i + 1 < argc) {
            recon_folder = argv[++i];
        } else if (arg == "--algorithm" && i + 1 < argc) {
            config.algorithm = argv[++i];
            std::transform(config.algorithm.begin(), config.algorithm.end(), 
                         config.algorithm.begin(), ::toupper);
        } else if (arg == "--iterations" && i + 1 < argc) {
            config.num_iterations = std::stoi(argv[++i]);
        } else if (arg == "--relaxation" && i + 1 < argc) {
            config.relaxation = std::stof(argv[++i]);
        } else if (arg == "--z-min" && i + 1 < argc) {
            config.recon_z_min = std::stoi(argv[++i]);
        } else if (arg == "--z-max" && i + 1 < argc) {
            config.recon_z_max = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            config.num_threads = std::stoi(argv[++i]);
        } else if (arg == "--no-log") {
            config.apply_log = false;
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Validate algorithm
    if (config.algorithm != "SART" && config.algorithm != "SIRT" && config.algorithm != "MART") {
        std::cerr << "Error: Unknown algorithm '" << config.algorithm << "'" << std::endl;
        std::cerr << "Valid options: SART, SIRT, MART" << std::endl;
        return 1;
    }
    
    // Set number of threads
    if (config.num_threads > 0) {
        omp_set_num_threads(config.num_threads);
    }
    
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Projection folder: " << proj_folder << std::endl;
    std::cout << "  Output folder: " << recon_folder << std::endl;
    std::cout << "  Algorithm: " << config.algorithm << std::endl;
    std::cout << "  Iterations: " << config.num_iterations << std::endl;
    std::cout << "  Relaxation: " << config.relaxation << std::endl;
    std::cout << "  Z range: [" << config.recon_z_min << ", " << config.recon_z_max << "]" << std::endl;
    std::cout << "  Threads: " << omp_get_max_threads() << std::endl;
    std::cout << "  Apply log: " << (config.apply_log ? "yes" : "no") << std::endl;
    std::cout << std::endl;
    
    // Check if input folder exists
    if (!fs::exists(proj_folder)) {
        std::cerr << "Error: Input folder not found: " << proj_folder << std::endl;
        std::cout << "\nCreating directory structure..." << std::endl;
        fs::create_directories(proj_folder);
        std::cout << "Please place your .tif projection files in: " << proj_folder << std::endl;
        return 1;
    }
    
    // Create output folder
    std::string output_folder = recon_folder + "/" + config.algorithm + "_reconstruction";
    fs::create_directories(output_folder);
    
    // Load projections
    std::vector<float> projections_data;
    if (!loadProjections(proj_folder, config, projections_data)) {
        std::cerr << "Error: Failed to load projections" << std::endl;
        return 1;
    }
    
    // Perform reconstruction
    if (!reconstructSliceBySlice(projections_data, config, output_folder)) {
        std::cerr << "Error: Reconstruction failed" << std::endl;
        return 1;
    }
    
    std::cout << "\n✓ All done!" << std::endl;
    
    return 0;
}