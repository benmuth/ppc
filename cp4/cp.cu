#include <cmath>
#include <cstdlib>
#include <cstring>
// #include <iostream>
// #include <ostream>
// #include <vector>

#include <cstdlib>
#include <cuda_runtime.h>
#include <iostream>

static inline void check(cudaError_t err, const char* context) {
  if (err != cudaSuccess) {
    std::cerr << "CUDA error: " << context << ": "
     << cudaGetErrorString(err) << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

#define CHECK(x) check(x, #x)

static inline int divup(int a, int b) { return (a + b - 1) / b; }

typedef float float8_t __attribute__((vector_size(8 * sizeof(float))));

// static inline int roundup(int a, int b) { return divup(a, b) * b; }

__global__ void mykernel(float *r, const float *d, int nx, int ny) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  if (i >= ny || j >= nx || i < j)
    return;

  float  v = 0;
  for (int k = 0; k < nx; ++k) {
    v += d[i * nx + k] * d[j * nx + k];
  }
  r[i + j * ny] = v;
}

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
  int data_size = nx * ny;

  float *normalized = (float *)malloc(data_size * sizeof(float));
  std::memset(normalized, 0, data_size * sizeof(float));

  // normalize input rows (center vector around origin)
  for (int j = 0; j < ny; ++j) {
    float row_sum = 0;
    int row_start = j * nx;
    for (int i = 0; i < nx; ++i) {
      row_sum += data[i + row_start];
    }
    float row_mean = row_sum / nx;

    for (int i = 0; i < nx; ++i) {
      normalized[i + row_start] = data[i + row_start] - row_mean;
    }

    float row_square_sum = 0;

    for (int i = 0; i < nx; ++i) {
      float val = normalized[i + row_start];
      row_square_sum += val * val;
    }

    float row_magnitude = sqrt(row_square_sum);

    for (int i = 0; i < nx; ++i) {
      normalized[i + row_start] = normalized[i + row_start] / row_magnitude;
    }
  }

  float *dGPU = NULL;
  float *rGPU = NULL;

  CHECK((cudaMalloc((void **)&dGPU, data_size * sizeof(float))));
  CHECK((cudaMalloc((void **)&rGPU, data_size * sizeof(float))));

  CHECK(cudaMemcpy(dGPU, normalized, data_size * sizeof(float),
                   cudaMemcpyHostToDevice));

  dim3 dimBlock(16, 16);
  dim3 dimGrid(divup(data_size, dimBlock.x), divup(data_size, dimBlock.y));
  mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, ny);
  CHECK((cudaGetLastError()));

  CHECK(cudaMemcpy(result, rGPU, data_size * sizeof(float),
                     cudaMemcpyDeviceToHost));

  CHECK(cudaFree(dGPU));
  CHECK(cudaFree(rGPU));

  cudaDeviceSynchronize();
  free(normalized);
}

