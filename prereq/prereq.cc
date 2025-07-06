struct Result {
  float avg[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- horizontal position: 0 <= x0 < x1 <= nx
- vertical position: 0 <= y0 < y1 <= ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
- output: avg[c]
*/
Result calculate(int ny, int nx, const float *data, int y0, int x0, int y1,
                 int x1) {

  Result result{{0.0f, 0.0f, 0.0f}};
  double r_sum = 0;
  double g_sum = 0;
  double b_sum = 0;
  
  for (int y = y0; y < y1; y++) {
    for (int x = x0; x < x1; x++) {
      int base_idx = 3 * x + 3 * nx * y;
      r_sum += data[base_idx + 0];
      g_sum += data[base_idx + 1];
      b_sum += data[base_idx + 2];
    }
  }
  int num_pixels = (x1 - x0) * (y1 - y0);
  result.avg[0] = r_sum / num_pixels;
  result.avg[1] = g_sum / num_pixels;
  result.avg[2] = b_sum / num_pixels;

  return result;
}
