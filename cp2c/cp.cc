#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ostream>
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

  int elems_per_padded_row =
      (nx + elems_per_vec - 1) / elems_per_vec * elems_per_vec;
  int vectors_per_padded_row = elems_per_padded_row / elems_per_vec;
  // int data_size = elems_per_padded_row * ny;

  double4_t zero_vec = {0.0, 0.0, 0.0, 0.0};
  std::vector<double4_t> vnormalized(vectors_per_padded_row * ny, zero_vec);
  std::vector<double4_t> vdata(vectors_per_padded_row * ny, zero_vec);

  // for (int j = 0; j < ny; ++j) {
  // }

  // normalize input rows (center vector of length 1 around origin)
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
  // correlation) first, iterate over pairs of rows (upper triangle)
  for (int i = 0; i < ny; ++i) {
    int row_i_start = i * vectors_per_padded_row;
    int i_ny = i * ny;
    for (int j = i; j < ny; ++j) {
      int row_j_start = j * vectors_per_padded_row;
      // iterate through the rows and sum the products of elements (dot product)
      double4_t correlation = {0.0, 0.0, 0.0, 0.0};

      // std::cerr << "row i start " << row_i_start << std::endl;
      // std::cerr << "row j start " << row_j_start << std::endl;
      for (int k = 0; k < vectors_per_padded_row; ++k) {
        correlation +=
            vnormalized[row_i_start + k] * vnormalized[row_j_start + k];
      }

      result[j + i_ny] =
          correlation[0] + correlation[1] + correlation[2] + correlation[3];
    }
  }
}
