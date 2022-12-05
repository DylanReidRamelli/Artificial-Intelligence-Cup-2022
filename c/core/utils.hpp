#ifndef UTILS_HPP
#define UTILS_HPP
#include <cmath>
#include <vector>

int compute_length(std::vector<int> solution,
                   std::vector<std::vector<int>> dist_matrix) {
  int total_length = 0;
  int starting_node = solution[0];
  int from_node = starting_node;
  for (int node : solution) {
    total_length += dist_matrix[from_node][node];
    from_node = node;
  }

  total_length += dist_matrix[starting_node][from_node];
  return total_length;
}

int distance_euc(std::vector<float> points_x, std::vector<float> points_y) {
  int rounding = 0;
  float x_i = points_x[0];
  float y_i = points_x[1];
  float x_j = points_y[0];
  float y_j = points_y[1];
  float distance = sqrt(pow((x_i - x_j), 2) + pow((y_i - y_j), 2));
  return floor(distance);
}

#endif