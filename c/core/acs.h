#ifndef ACS_H
#define ACS_H
#include "io_tsp.h"

typedef struct {
  int start_node;
  int index;
  int *trail;
  int *visited;
  // int trail[130];
  double cost;
} ANT;

void nearest_neighbor(int *sol, int distances[dim_][dim_], int *new_sol) {
  int start_node = rand() % dim_; // index of randomly selected node
  // int start_node = find_index(sol,0);
  int visited[dim_];
  for (int i = 0; i < dim_; i++) {
    visited[i] = 0;
  }
  visited[sol[start_node]] = 1;
  new_sol[0] = sol[start_node];
  int nearest;
  int value;
  int index = sol[start_node];
  for (int j = 1; j < dim_; j++) {
    value = distances[0][0];
    for (size_t k = 0; k < dim_; k++) {
      if (distances[index][k] <= value && visited[k] == 0) {
        nearest = k;
        value = distances[index][nearest];
      }
    }
    // printf("%d\n", index);
    index = nearest;
    visited[index] = 1;
    new_sol[j] = index;
  }
}

int check_valid(int *sol) {
  int visited[dim_];
  for (size_t i = 0; i < dim_; i++) {
    visited[i] = 0;
  }
  for (size_t j = 0; j < dim_; j++) {
    int node = sol[j];
    visited[node] = visited[node] += 1;
  }
  for (int k = 0; k < dim_; k++) {
    if (visited[k] != 1) {
      // printf("%d\t%d\t%d\n",sol[0],k,visited[k] );
      // printf("%s\n", "Invalid Tour");
      return 1;
    }
  }
  return 0;
}
void global_pheromone_update(ANT *best, double pheromone[dim_][dim_],
                             double alpha) {
  int node1, node2;
  double phi, delta;
  delta = 1.0 / best->cost;
  for (size_t i = 0; i < dim_ - 1; i++) {
    node1 = best->trail[i];
    node2 = best->trail[i + 1];
    phi = pheromone[node1][node2];
    pheromone[node1][node2] = ((1 - alpha) * phi) + (alpha * delta);
    pheromone[node2][node1] = pheromone[node1][node2];
  }
  node1 = best->trail[0];
  node2 = best->trail[dim_ - 1];
  phi = pheromone[node1][node2];
  pheromone[node1][node2] = ((1 - alpha) * phi) + (alpha * delta);
  pheromone[node2][node1] = pheromone[node1][node2];
}

void acs(ANT **nest, int swarm_count, int distances[dim_][dim_],
         double pheromone[dim_][dim_], int *shortest_path) {
  // declare variables
  ANT best;
  int i, j; // loop indices
  int node, q, valid, pos;
  double tau, min, alpha, t;
  alpha = 0.1;
  t = 0.1; // evaporation
  min = 99999999;
  ANT *a;
  int flag = 0; // down
  /*
   * Initialise ants.
   * Put m ants on randomly selected cities
   */
  clean_ants(nest, swarm_count);
  // make tour for each ant
  for (i = 0; i < swarm_count; i++) { // for the current ant
    flag = 0;
    a = nest[i];
    pos = a->start_node;
    printf("%d\n", pos);
    for (j = 0; j < dim_; j++) { // pick cities
      q = rand() % 100;
      if (q <= 94) {
        node = exploit(a, pos, a->visited, pheromone, distances);
        // printf("%d\t", node);

      } else {
        node = explore(a, pos, pheromone, distances);
      }
      // add the node
      if (a->visited[node] == 0 && node != -1 && a->index < dim_) {
        a->trail[a->index] = node;
        tau = pheromone[pos][node];
        // local pheromone update (symmetric)
        pheromone[pos][node] = (1 - alpha) * tau + (alpha * t);
        pheromone[node][pos] = (1 - alpha) * tau + (alpha * t);
        a->index++;
        a->visited[node] = 1;
        pos = node;
        node = -1;
      } else {
        flag = 1;
      }
    }
    // update ant details
    valid = check_valid(a->trail);
    if (valid == 0 && flag != 1) {
      a->cost = calculate_solution_cost(a->trail, dim_, distances);
    } else {
      a->cost = 99999999;
    }
  }
  for (i = 0; i < swarm_count; i++) {
    if (nest[i]->cost < min) {
      min = nest[i]->cost;
      // best = *nest[i];
      for (size_t j = 0; j < dim_; j++) {
        best.trail[j] = nest[i]->trail[j];
        shortest_path[j] = nest[i]->trail[j];
      }
    }
  }
  global_pheromone_update(&best, pheromone, alpha);
}

int exploit(ANT *a, int pos, int *visited, double pheromone[dim_][dim_],
            int distances[dim_][dim_]) {
  double tau, eta, value;
  // printf("%d\n", a->index);
  // int pos = a->trail[a->index-1];
  double beta = 2.0;
  double max = -1;
  int node = -1;
  for (int k = 0; k < dim_; k++) {
    tau = pheromone[pos][k];
    eta = pow((1.0 / (double)distances[pos][k]), beta);
    value = tau * eta;
    if (value > max && visited[k] == 0 && pos != k) {
      max = value;
      node = k;
    }
  }
  return node;
}
int explore(ANT *a, int pos, double pheromone[dim_][dim_],
            int distances[dim_][dim_]) {
  double tau, eta, value, rnd, d;
  double sigma = 0;
  double S[dim_];
  double beta = 2.0;
  int node = -1;
  // int pos = a->trail[a->index-1];
  // calculate probability
  for (int i = 0; i < dim_; i++) {
    if (a->visited[i] == 0) {
      tau = pheromone[pos][i];
      eta = pow((1.0 / (double)distances[pos][i]), beta);
      value = tau * eta;
      S[i] = value;
      sigma += value;
    } else {
      S[i] = 0;
    }
  }
  for (int i = 0; i < dim_; i++) {
    S[i] = S[i] / sigma;
  }

  // Pick a random edge
  rnd = (double)rand() / ((double)RAND_MAX);
  // printf("%lf\n", rnd);
  d = 0;
  for (int i = 0; i < dim_; i++) {
    d += S[i];
    if (rnd < d && S[i] != 0) {
      //    printf("%lf\t %lf\t %d\n", rnd,d,i);
      node = i;
      return node;
    }
  }
  return node;
}
void clean_ants(ANT **nest, int swarm_count) {
  ANT *a1;
  for (int i = 0; i < swarm_count; i++) {
    a1 = nest[i];
    a1->start_node = rand() % dim_;
    a1->index = 0;
    a1->cost = 9999;
    for (int j = 0; j < dim_; j++) {
      a1->visited[j] = 0;
      a1->trail[j] = 0;
    }
    a1->trail[a1->index] = a1->start_node;
    a1->visited[a1->start_node] = 1;
    a1->index++;
  }
}

#endif