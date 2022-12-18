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

int distance_euc(std::vector<int> points_x, std::vector<int> points_y) {
  int rounding = 0;
  int x_i = points_x[0];
  int y_i = points_x[1];
  int x_j = points_y[0];
  int y_j = points_y[1];
  float distance = sqrt(pow((x_i - x_j), 2) + pow((y_i - y_j), 2));
  return floor(distance);
}

void calculate_distance(std::vector<std::vector<int>> nodes,
                        std::vector<std::vector<int>> distances, int dim) {
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < dim; j++) {
      if (j == i) {
        distances[i][j] = 9999;
        distances[j][i] = 9999;
      } else {
        double x_sq = pow(nodes[0][i] - nodes[0][j], 2);
        double y_sq = pow(nodes[1][i] - nodes[1][j], 2);
        distances[i][j] = sqrt(x_sq + y_sq);
        distances[j][i] = sqrt(x_sq + y_sq);
      }
    }
  }
}

int calculate_solution_cost(std::vector<int> shortest_path,
                            std::vector<std::vector<int>> distances, int dim) {
  int sum = 0;
  for (size_t i = 0; i < dim - 1; i++) {
    sum += distances[shortest_path[i]][shortest_path[i + 1]];
  }
  sum += distances[shortest_path[dim - 1]][shortest_path[0]];
  return sum;
}

#endif