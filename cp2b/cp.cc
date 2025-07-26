#include <cmath>
#include <cstdlib>

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

  double *normalized = (double *)malloc(data_size * sizeof(double));
  #pragma omp parallel for
  for (int i = 0; i < data_size; ++i) {
    normalized[i] = 0;
  }

  // normalize input rows (center vector around origin)

  #pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    double row_sum = 0;
    for (int i = 0; i < nx; ++i) {
      row_sum += data[i + j * nx];
    }
    double row_mean = row_sum / nx;

    for (int i = 0; i < nx; ++i) {
      normalized[i + (j * nx)] = data[i + (j * nx)] - row_mean;
    }
  }

  // normalize input rows (scale vector to length of 1)
  #pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    double row_square_sum = 0;
    for (int i = 0; i < nx; ++i) {
      double val = normalized[i + j * nx];
      row_square_sum += val * val;
    }
    double row_magnitude = sqrt(row_square_sum);

    for (int i = 0; i < nx; ++i) {
      normalized[i + (j * nx)] = normalized[i + (j * nx)] / row_magnitude;
    }
  }

  // calculate upper triangle of matrix product (normalized matrix) x
  // (transpose)
  // this is taking the dot product of each pairwise row (which gives the
  // correlation) first, iterate over pairs of rows (upper triangle)
  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < ny; ++i) {
    for (int j = i; j < ny; ++j) {
      double correlation = 0;
      // second, iterate through rows and sum products of elements (dot product)
      for (int k = 0; k < nx; ++k) {
        correlation += normalized[i * nx + k] * normalized[j * nx + k];
      }
      result[j + i * ny] = correlation;
    }
  }

  free(normalized);
}
