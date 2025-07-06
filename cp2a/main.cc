#include <iostream>
#include <random>
#include <string>
#include <chrono>
#include <cstdlib>

// Forward declaration of the correlate function from cp.cc
void correlate(int ny, int nx, const float *data, float *result);

int main() {
    // Hard-coded benchmark parameters from benchmarks/2.txt
    int ny = 4000, nx = 1000, seed = 5;
    
    std::cout << "Running benchmark with ny=" << ny << ", nx=" << nx << ", seed=" << seed << std::endl;
    
    // Generate random data
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    
    float *data = new float[ny * nx];
    for (int i = 0; i < ny * nx; ++i) {
        data[i] = dis(gen);
    }
    
    // Allocate result matrix (upper triangle only)
    float *result = new float[ny * ny];
    
    // Time the correlation computation
    auto start = std::chrono::high_resolution_clock::now();
    correlate(ny, nx, data, result);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Correlation computation took " << duration.count() << " ms" << std::endl;
    
    // Print a few sample results for verification
    std::cout << "Sample correlations:" << std::endl;
    for (int i = 0; i < std::min(5, ny); ++i) {
        for (int j = i; j < std::min(i + 5, ny); ++j) {
            std::cout << "corr(" << i << "," << j << ") = " << result[j + i * ny] << " ";
        }
        std::cout << std::endl;
    }
    
    delete[] data;
    delete[] result;
    
    return 0;
}