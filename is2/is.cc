// #include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>

// precomputation
// compute color sums for all rectangles anchored at (0,0) using a single-pass
// summed area table algorithm essentially a 2d prefix sum

// precompute color squared sum ( (a_i)^2 )

// iterate through all combinations of rectangles and keep track of minimum SSE
// SSE = n * p - 2(a_i)p + (a_i)^2

struct precomp_values {
  double *areas;
  double *color_squared_sum;
};

size_t index(int c, int x, int y, int nx) { return c + 3 * x + 3 * nx * y; }

precomp_values make_summed_area(int ny, int nx, const float *data) {
  double *areas = (double *)calloc(nx * ny * 3, sizeof(double));
  double *color_squared_sum = (double *)calloc(nx * ny * 3, sizeof(double));

  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      for (int c = 0; c < 3; c++) {
        int i = index(c, x, y, nx);

        areas[i] = data[index(c, x, y, nx)];

        if (x - 1 >= 0) {
          areas[i] += areas[index(c, x - 1, y, nx)];
        }

        if (y - 1 >= 0) {
          areas[i] += areas[index(c, x, y - 1, nx)];
        }

        if (x - 1 >= 0 && y - 1 >= 0) {
          areas[i] -= areas[index(c, x - 1, y - 1, nx)];
        }
      }
    }
  }

  return precomp_values{.areas = areas, .color_squared_sum = color_squared_sum};
}

struct Result {
  // upper left corner
  int y0;
  int x0;
  // lower right corner (exclusive)
  int y1;
  int x1;
  // color of background
  float outer[3];
  // color of rectangle
  float inner[3];
};

struct Rectangle {
  int x0;
  int y0;
  int x1;
  int y1;
  float color[3];
};

double get_summed_area(double *areas, int x, int y, int c, int width,
                       int height, int nx) {
  double sum = areas[index(c, x + width, y + height, nx)];
  sum -= areas[index(c, x + width, y, nx)];
  sum -= areas[index(c, x, y + height, nx)];
  sum += areas[index(c, x, y, nx)];
  return sum;
}

float error(const float *data, int nx, int ny, Result r);

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/
Result segment(int ny, int nx, const float *data) {
  Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

  // precomputation
  Rectangle *rects = (Rectangle *)malloc(nx * ny * sizeof(Rectangle));
  for (int y = 0; y < ny; ++y) {
    float color[3] = {0.0, 0.0, 0.0};
    for (int x = 0; x < nx; ++x) {
      // indices of the rectangle and the rectangle above
      int i = (y * nx) + x;
      int prev_row_i = ((y - 1) * nx) + x;

      // calculate the 2D prefix sum of color for each rectangle (easy to find
      // average color for a rectangle later)

      // 00 01 02
      // 10 11 12

      // sum(11) = 11 + 10 + 01 - 00
      if (prev_row_i >= 0 && i > 0) {
        for (int c = 0; c < 3; c++) {
          color[c] +=
              rects[prev_row_i].color[c] + data[c + (3 * x) + ((3 * nx) * y)];
        }
      } else {
        for (int c = 0; c < 3; ++c) {
          color[c] += data[c + (3 * x) + ((3 * nx) * y)];
        }
      }
      // float prev_color[3];
      // prev_color[0] = rects[i - 1].color[0];
      // prev_color[1] = rects[i - 1].color[1];
      // prev_color[2] = rects[i - 1].color[2];
      // std::cout << color[0] << std::endl;
      rects[i] = Rectangle{.x0 = 0,
                           .y0 = 0,
                           .x1 = x,
                           .y1 = y,
                           .color = {color[0], color[1], color[2]}};
    }
  }

  // // print precomputed rectangles
  // std::cout << std::fixed << std::setprecision(2);
  // for (int i = 0; i < nx * ny; ++i) {
  //   std::cout << "[" << rects[i].x0 << "," << rects[i].y0 << "," <<
  //   rects[i].x1
  //             << "," << rects[i].y1 << "] ";

  //   for (int c = 0; c < 3; ++c) {
  //     std::cout << rects[i].color[c] << " ";
  //   }
  //   std::cout << " | ";

  //   if ((i + 1) % nx == 0) {
  //     std::cout << std::endl;
  //   }
  // }

  float min_error = 1e10;
  // find best rectangle
  // try every rectangle?
  for (int y0 = 0; y0 < ny; ++y0) {
    for (int x0 = 0; x0 < nx; ++x0) {
      for (int y1 = y0; y1 < ny; ++y1) {
        for (int x1 = x0; x1 < nx; ++x1) {
          float err = 0.0;

          Result candidate = get_summed_area(precomp_values->areas, int x, int y, int c,
                                             int width, int height, int nx);
          err += error(data, nx, ny, &candidate);
          if (err < min_error) {
            result = candidate;
            min_error = err;
          }
        }
      }
    }
  }

  std::cout << std::endl;

  return result;
}

float error(const float *data, int nx, int ny, Result *r) {
  // error
  // float outside_color[3] = {150, 125, 112};
  // float inside_color[3] = {112, 125, 150};
  // Rectangle rect = {
  //     .x0 = 0, .y0 = 0, .x1 = 25, .y1 = 25, .color = {100, 100, 100}};

  float sum_squared_error = 0.0;

  // SSE = n * p - 2(a_i)p + (a_i)^2

  // for (int x = 0; x < nx; ++x) {
  //   for (int y = 0; y < ny; ++y) {
  //     for (int c = 0; c < 3; ++c) {
  //       float error = 0.0;
  //       float point_color = data[c + 3 * x + 3 * nx * y];

  //       // outside
  //       if (x < r->x0 || x > (r->x1 - 1) || y < r->y0 || (y > r->y1 - 1)) {
  //         error = r->outer[c] - point_color;
  //       } else { // inside
  //         error = r->inner[c] - point_color;
  //       }
  //       sum_squared_error += error * error;
  //     }
  //   }
  // }

  return sum_squared_error;
}

// Rectangle createRect(Rectangle *precomputed_rects, int size, Rectangle
// want_rect) {
//   Rectangle upper_left;
//   Rectangle lower_left;
//   Rectangle upper_right;
//   Rectangle lower_right;

//   for (int i = 0; i < size; ++i) {

//   }
// }
