#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));

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
  constexpr int stride = 4;

  int elems_per_padded_row =
      (nx + elems_per_vec - 1) / elems_per_vec * elems_per_vec;
  int vectors_per_padded_row = elems_per_padded_row / elems_per_vec;

  double4_t zero_vec = {0.0, 0.0, 0.0, 0.0};
  std::vector<double4_t> vnormalized(vectors_per_padded_row * ny, zero_vec);
  std::vector<double4_t> vdata(vectors_per_padded_row * ny, zero_vec);

  // normalize input rows (center vector of length 1 around origin)
  #pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    double4_t sum_vec = {0.0, 0.0, 0.0, 0.0};
    int row_start_vec = j * vectors_per_padded_row;
    // pack data vectors
    for (int i = 0; i < nx; ++i) {
      int vector_idx = row_start_vec + i / elems_per_vec;
      int elem_idx = i % elems_per_vec;
      vdata[vector_idx][elem_idx] = data[i + j * nx];
    }

    // sum vectors
    for (int i = 0; i < vectors_per_padded_row; ++i) {
      sum_vec += vdata[i + row_start_vec];
    }
    double row_sum = sum_vec[0] + sum_vec[1] + sum_vec[2] + sum_vec[3];
    // calculate mean
    double row_mean = row_sum / nx;

    // subtract mean
    double4_t mean_vec = {row_mean, row_mean, row_mean, row_mean};
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
    double4_t row_square_sum_vec = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < vectors_per_padded_row; ++i) {
      double4_t val = vnormalized[i + row_start_vec];
      row_square_sum_vec += val * val;
    }

    double row_square_sum = row_square_sum_vec[0] + row_square_sum_vec[1] +
                            row_square_sum_vec[2] + row_square_sum_vec[3];
    double row_magnitude = sqrt(row_square_sum);

    double4_t row_magnitude_vec = {1 / row_magnitude, 1 / row_magnitude,
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
      double4_t correlation0 = {0.0, 0.0, 0.0, 0.0};
      double4_t correlation1 = {0.0, 0.0, 0.0, 0.0};
      double4_t correlation2 = {0.0, 0.0, 0.0, 0.0};
      double4_t correlation3 = {0.0, 0.0, 0.0, 0.0};

      int row_j0_start = j * vectors_per_padded_row;
      int row_j1_start = (j + 1 < ny) ? (j + 1) * vectors_per_padded_row : 0;
      int row_j2_start = (j + 2 < ny) ? (j + 2) * vectors_per_padded_row : 0;
      int row_j3_start = (j + 3 < ny) ? (j + 3) * vectors_per_padded_row : 0;

      // Compute dot products for up to stride rows simultaneously
      for (int k = 0; k < vectors_per_padded_row; ++k) {
        double4_t row_i_vec = vnormalized[row_i_start + k];
        correlation0 += row_i_vec * vnormalized[row_j0_start + k];
        if (j + 1 < ny) correlation1 += row_i_vec * vnormalized[row_j1_start + k];
        if (j + 2 < ny) correlation2 += row_i_vec * vnormalized[row_j2_start + k];
        if (j + 3 < ny) correlation3 += row_i_vec * vnormalized[row_j3_start + k];
      }

      // Store results
      result[j + i_ny] = correlation0[0] + correlation0[1] + correlation0[2] + correlation0[3];
      if (j + 1 < ny) result[(j + 1) + i_ny] = correlation1[0] + correlation1[1] + correlation1[2] + correlation1[3];
      if (j + 2 < ny) result[(j + 2) + i_ny] = correlation2[0] + correlation2[1] + correlation2[2] + correlation2[3];
      if (j + 3 < ny) result[(j + 3) + i_ny] = correlation3[0] + correlation3[1] + correlation3[2] + correlation3[3];
    }
  }
}
