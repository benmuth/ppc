#include <cmath>
#include <cstdlib>
#include <cstring>

// operations
// subtract mean from each row
//

// x = [0, 1, 2, 3]
//

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
  int stride_length = 16;

  int data_size = nx * ny;

  double *normalized = (double *)malloc(data_size * sizeof(double));
  std::memset(normalized, 0, data_size * sizeof(double));

  // normalize input rows (center vector around origin)
  for (int j = 0; j < ny; ++j) {
    double row_sum = 0;
    int row_start = j * nx;
    for (int i = 0; i < nx; ++i) {
      row_sum += data[i + row_start];
    }
    double row_mean = row_sum / nx;

    for (int i = 0; i < nx; ++i) {
      normalized[i + row_start] = data[i + row_start] - row_mean;
    }

    double row_square_sum = 0;

    for (int i = 0; i < nx; ++i) {
      double val = normalized[i + row_start];
      row_square_sum += val * val;
    }

    double row_magnitude = sqrt(row_square_sum);

    for (int i = 0; i < nx; ++i) {
      normalized[i + row_start] = normalized[i + row_start] / row_magnitude;
    }
  }

  // calculate upper triangle of matrix product (normalized matrix) x
  // (transpose)
  // this is taking the dot product of each pairwise row (which gives the
  // correlation) first, iterate over pairs of rows (upper triangle)
  for (int i = 0; i < ny; ++i) {
    int row_i_start = i * nx;
    int i_ny = i * ny;
    for (int j = i; j < ny; ++j) {
      int row_j_start = j * nx;
      // second, iterate through the rows and sum the products of elements (dot
      // product)
      double correlation = 0;
      double correlation1 = 0;
      double correlation2 = 0;
      double correlation3 = 0;
      double correlation4 = 0;
      double correlation5 = 0;
      double correlation6 = 0;
      double correlation7 = 0;
      double correlation8= 0;
      double correlation9 = 0;
      double correlation10 = 0;
      double correlation11 = 0;
      double correlation12 = 0;
      double correlation13 = 0;
      double correlation14 = 0;
      double correlation15 = 0;
      for (int k = 0; k < nx / stride_length * stride_length;
           k += stride_length) {
        correlation +=
            normalized[row_i_start + k] * normalized[row_j_start + k];
        correlation1 +=
            normalized[row_i_start + k + 1] * normalized[row_j_start + k + 1];
        correlation2 +=
            normalized[row_i_start + k + 2] * normalized[row_j_start + k + 2];
        correlation3 +=
            normalized[row_i_start + k + 3] * normalized[row_j_start + k + 3];
        correlation4 +=
            normalized[row_i_start + k + 4] * normalized[row_j_start + k + 4];
        correlation5 +=
            normalized[row_i_start + k + 5] * normalized[row_j_start + k + 5];
        correlation6 +=
            normalized[row_i_start + k + 6] * normalized[row_j_start + k + 6];
        correlation7 +=
            normalized[row_i_start + k + 7] * normalized[row_j_start + k + 7];
        correlation8 +=
            normalized[row_i_start + k + 8] * normalized[row_j_start + k + 8];
        correlation9  +=
            normalized[row_i_start + k + 9] * normalized[row_j_start + k + 9];
        correlation10  +=
            normalized[row_i_start + k + 10] * normalized[row_j_start + k + 10];
        correlation11  +=
            normalized[row_i_start + k + 11] * normalized[row_j_start + k + 11];
        correlation12  +=
            normalized[row_i_start + k + 12] * normalized[row_j_start + k + 12];
        correlation13  +=
            normalized[row_i_start + k + 13] * normalized[row_j_start + k + 13];
        correlation14  +=
            normalized[row_i_start + k + 14] * normalized[row_j_start + k + 14];
        correlation15  +=
            normalized[row_i_start + k + 15] * normalized[row_j_start + k + 15];
      }

      for (int k = nx / stride_length * stride_length; k < nx; k++) {
        correlation +=
            normalized[row_i_start + k] * normalized[row_j_start + k];
      }

      result[j + i_ny] = correlation + correlation1 + correlation2 +
                         correlation3 + correlation4 + correlation5 +
                         correlation6 + correlation7 + correlation8+ correlation9 + correlation10 +
                         correlation11 + correlation12 + correlation13 +
                         correlation14 + correlation15;
    }
  }

  free(normalized);
}
