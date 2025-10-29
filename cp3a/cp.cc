#include <cmath>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>
#include <vector>

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
  constexpr int elems_per_vec = 4;
  constexpr int stride = 2;

  int elems_per_padded_row =
      ((nx + elems_per_vec - 1) / elems_per_vec) * elems_per_vec;

  int vectors_per_padded_row = elems_per_padded_row / elems_per_vec;

  __m256d zero_vec = _mm256_setzero_pd();
  std::vector<__m256d> vnormalized(vectors_per_padded_row * ny, zero_vec);
  std::vector<__m256d> vdata(vectors_per_padded_row * ny, zero_vec);

// normalize input rows (center vector of length 1 around origin)
#pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    __m256d sum_vec = _mm256_setzero_pd();
    int row_start_vec = j * vectors_per_padded_row;
    // pack data vectors
    for (int i = 0; i < nx; i += elems_per_vec) {
      int vector_idx = row_start_vec + i / elems_per_vec;
      vdata[vector_idx] =
          _mm256_set_pd(data[(i) + j * nx], data[(i + 1) + j * nx],
                        data[(i + 2) + j * nx], data[(i + 3) + j * nx]);
    }

    // sum vectors
    for (int i = 0; i < vectors_per_padded_row; ++i) {
      _mm256_add_pd(sum_vec, vdata[i + row_start_vec]);
    }
    double row_sum = sum_vec[0] + sum_vec[1] + sum_vec[2] + sum_vec[3];
    // calculate mean
    double row_mean = row_sum / nx;

    // subtract mean
    __m256d mean_vec = {row_mean, row_mean, row_mean, row_mean};
    for (int i = 0; i < vectors_per_padded_row; ++i) {
      vnormalized[i + row_start_vec] = vdata[i + row_start_vec] - mean_vec;
    }

    // set padded elements to 0
    if (nx % elems_per_vec != 0) {
      int padding = nx % elems_per_vec;
      for (int k = padding; k < elems_per_vec; ++k) {
        vnormalized[row_start_vec + vectors_per_padded_row - 1][k] = 0.0;
      }
    }

    // square elems
    __m256d row_square_sum_vec = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < vectors_per_padded_row; ++i) {
      __m256d val = vnormalized[i + row_start_vec];
      row_square_sum_vec += val * val;
    }

    double row_square_sum = row_square_sum_vec[0] + row_square_sum_vec[1] +
                            row_square_sum_vec[2] + row_square_sum_vec[3];
    double row_magnitude = sqrt(row_square_sum);

    __m256d row_magnitude_vec = {1 / row_magnitude, 1 / row_magnitude,
                                 1 / row_magnitude, 1 / row_magnitude};
    // normalize
    for (int i = 0; i < vectors_per_padded_row; ++i) {
      vnormalized[i + row_start_vec] *= row_magnitude_vec;
    }
  }

// calculate upper triangle of matrix product (normalized matrix) x
// (transpose)
// this is taking the dot product of each pairwise row (which gives the
// correlation)
// first, iterate over pairs of rows (upper triangle)
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < ny; ++i) {
    int row_i_start = i * vectors_per_padded_row;
    int i_ny = i * ny;

    // Process stride rows at a time
    for (int j = i; j < ny; j += stride) {
      __m256d correlation0 = {0.0, 0.0, 0.0, 0.0};
      __m256d correlation1 = {0.0, 0.0, 0.0, 0.0};

      int row_j0_start = j * vectors_per_padded_row;
      int row_j1_start = (j + 1 < ny) ? (j + 1) * vectors_per_padded_row : 0;

      // Compute dot products for up to stride rows simultaneously
      for (int k = 0; k < vectors_per_padded_row; ++k) {
        __m256d row_i_vec = vnormalized[row_i_start + k];
        correlation0 += row_i_vec * vnormalized[row_j0_start + k];
        if (j + 1 < ny)
          correlation1 += row_i_vec * vnormalized[row_j1_start + k];
      }

      // Store results
      result[j + i_ny] =
          correlation0[0] + correlation0[1] + correlation0[2] + correlation0[3];
      if (j + 1 < ny)
        result[(j + 1) + i_ny] = correlation1[0] + correlation1[1] +
                                 correlation1[2] + correlation1[3];
    }
  }
}
