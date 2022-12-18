#ifndef UTILS_H
#define UTILS_H

void calculate_distance(double nodes[2][dim], int distances[dim][dim]) {
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
int calculate_solution_cost(int *sol, int size, int distances[dim][dim]) {
  int sum = 0;
  for (size_t i = 0; i < dim - 1; i++) {
    sum += distances[sol[i]][sol[i + 1]];
  }
  sum += distances[sol[dim - 1]][sol[0]];
  return sum;
}

void swap_places(int *sol, int i, int j) {
  int temp = sol[j];
  sol[j] = sol[i];
  sol[i] = temp;
}
int find_index(int *sol, int val) {
  for (size_t i = 0; i < dim; i++) {
    if (sol[i] == val) {
      return i;
    }
  }
  return -1;
}

#endif