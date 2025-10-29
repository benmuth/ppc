#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>

struct Result {
  int y0;
  int x0;
  int y1;
  int x1;
  float outer[3];
  float inner[3];
};

struct Rectangle {
  int x0;
  int y0;
  int x1;
  int y1;
  float color[3];
};

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
      int i = (y * nx) + x;
      int prev_row_i = ((y - 1) * nx) + x;
      if (prev_row_i >= 0 && i > 0) {
        for (int c = 0; c < 3; ++c) {
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

  // print precomputed rectangles
  std::cout << std::fixed << std::setprecision(2);
  for (int i = 0; i < nx * ny; ++i) {
    std::cout << "[" << rects[i].x0 << "," << rects[i].y0 << "," << rects[i].x1
              << "," << rects[i].y1 << "] ";

    for (int c = 0; c < 3; ++c) {
      std::cout << rects[i].color[c] << " ";
    }
    std::cout << " | ";

    if ((i + 1) % nx == 0) {
      std::cout << std::endl;
    }
  }

  float min_error = 1e10;
  // find best rectangle
  // try every rectangle?
  for (int y0 = 0; y0 < ny; ++y0) {
    for (int x0 = 0; x0 < nx; ++x0) {
      for (int y1 = 0; y1 < ny; ++y1) {
        for (int x1 = 0; x1 < nx; ++x1) {
          float err = 0.0;
          Result candidate = {};
          err += error(data, nx, ny, candidate);
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

float error(const float *data, int nx, int ny, Result r) {
  // error
  // float outside_color[3] = {150, 125, 112};
  // float inside_color[3] = {112, 125, 150};
  // Rectangle rect = {
  //     .x0 = 0, .y0 = 0, .x1 = 25, .y1 = 25, .color = {100, 100, 100}};

  float sum_squared_error = 0.0;

  for (int x = 0; x < nx; ++x) {
    for (int y = 0; y < ny; ++y) {
      for (int c = 0; c < 3; ++c) {
        float error = 0.0;
        float point_color = data[c + 3 * x + 3 * nx * y];

        // outside
        if (x < r.x0 || x > (r.x1 - 1) || y < r.y0 || (y > r.y1 - 1)) {
          error = r.outer[c] - point_color;
        } else { // inside
          error = r.inner[c] - point_color;
        }
        sum_squared_error += error * error;
      }
    }
  }

  return sum_squared_error;
}

Rectangle createRect(Rectangle *precomputed_rects, int size, Rectangle want_rect) {
  Rectangle upper_left;
  Rectangle lower_left;
  Rectangle upper_right;
  Rectangle lower_right;

  for (int i = 0; i < size; ++i) {
    
  }
}
