#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int dim;
// static int best;
static int *distances;
static double *inv_dist;
static double *pheromone;
static double *nodes;
static double t;
static int *tour;
static int tour_cost;

typedef struct {
  int start_node;
  int index;
  int *trail;
  int *visited;
  double cost;
} ANT;

// function declarations
int main(int argc, char const *argv[]);
FILE *read_words(FILE *f);
FILE *read_coords(FILE *f);
void calculate_distance();
int calculate_solution_cost(int *sol, int size);
void swap_places(int *sol, int i, int j);
void nearest_neighbor(int *sol, int *new_sol);
int check_valid(int *sol);
void acs(int swarm_count, int *shortest_path);
void global_pheromone_update(double alpha);
int *two_opt(int *sol);
int *two_opt_first(int *sol, int *distances);
void two_opt_swap(int *s, int i, int k);
int *hill_climbing(int *sol, int *distances);
int compute_gain(int *distances, int *sol, int i, int j);
void exchange(int *sol, int best_i, int best_j);

FILE *read_words(FILE *f) {
  char string[100];
  fscanf(f, "%s", string);
  while (strcmp(string, "NODE_COORD_SECTION") != 0) {
    if (strcmp(string, "DIMENSION:") == 0 || strcmp(string, "DIMENSION") == 0) {
      fscanf(f, "%s", string);
      if (strcmp(string, ":") == 0) {
        fscanf(f, "%s", string);
      }
      dim = atoi(string);
    }
    // if (strcmp(string,"BEST_KNOWN:")==0 || strcmp(string,"BEST_KNOWN")==0 ) {
    //     fscanf(f, "%s",string);
    //     if (strcmp(string,":")==0 ) {
    //         fscanf(f, "%s",string);
    //     }
    //     best = atoi(string);
    // }
    fscanf(f, "%s", string);
  }
  return f;
}
FILE *read_coords(FILE *f) {
  int index[1];
  double x[1];
  double y[1];
  for (size_t i = 0; i < dim; i++) {
    fscanf(f, "%d", index);
    fscanf(f, "%lf", x);
    fscanf(f, "%lf", y);
    *(nodes + i) = x[0];
    *(nodes + dim + i) = y[0];
  }
  return f;
}

void calculate_distance() {
  double x_sq, y_sq, value;

  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < dim; j++) {
      if (j == i) {
        *(distances + (i * dim + j)) = INT_MAX;
        *(distances + (j * dim + i)) = INT_MAX; // isnt this the same spot
      } else {
        value = *(nodes + i) - *(nodes + j);
        x_sq = value * value;
        // x_sq = pow(*(nodes+i) - *(nodes+j), 2);
        value = *(nodes + dim + i) - *(nodes + dim + j);
        y_sq = value * value;
        // y_sq = pow(*(nodes+dim+i) - *(nodes+dim+j), 2);
        value = round(sqrt(x_sq + y_sq));
        *(distances + (i * dim + j)) = value;
        *(distances + (j * dim + i)) = value;
      }
      *(inv_dist + (i * dim + j)) = 1.0 / *(distances + (i * dim + j));
      *(inv_dist + (j * dim + i)) = 1.0 / *(distances + (j * dim + i));
    }
  }
}
int calculate_solution_cost(int *sol, int size) {
  int sum = 0;
  for (size_t i = 0; i < dim - 1; i++) {
    sum += *(distances + (sol[i] * dim + sol[i + 1]));
  }
  sum += *(distances + (sol[dim - 1] * dim + sol[0]));
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

void nearest_neighbor(int *sol, int *new_sol) {
  int start_node = rand() % dim; // index of randomly selected node
  // int start_node = find_index(sol,0);
  int visited[dim];
  for (int i = 0; i < dim; i++) {
    visited[i] = 0;
  }
  visited[sol[start_node]] = 1;
  new_sol[0] = sol[start_node];
  int nearest;
  int value;
  int index = sol[start_node];
  for (int j = 1; j < dim; j++) {
    value = *distances;
    for (size_t k = 0; k < dim; k++) {
      if (*(distances + (index * dim + k)) <= value && visited[k] == 0) {
        nearest = k;
        value = *(distances + (index * dim + nearest));
      }
    }
    index = nearest;
    visited[index] = 1;
    new_sol[j] = index;
  }
}

int check_valid(int *sol) {
  int visited[dim];
  for (size_t i = 0; i < dim; i++) {
    visited[i] = 0;
  }
  for (size_t j = 0; j < dim; j++) {
    int node = sol[j];
    visited[node] = visited[node] += 1;
  }
  for (int k = 0; k < dim; k++) {
    if (visited[k] != 1) {
      return 1;
    }
  }
  return 0;
}
void global_pheromone_update(double alpha) {
  int node1, node2;
  double phi, delta;
  delta = 1.0 / tour_cost;
  for (size_t i = 0; i < dim - 1; i++) {
    node1 = tour[i];
    node2 = tour[i + 1] % dim;
    phi = *(pheromone + (node1 * dim + node2));
    *(pheromone + (node1 * dim + node2)) =
        ((1 - alpha) * phi) + (alpha * delta);
    *(pheromone + (node2 * dim + node1)) = *(pheromone + (node1 * dim + node2));
  }
  // node1 = tour[0];
  // node2 = tour[dim-1];
  // phi = *(pheromone+ (node1*dim + node2));
  // *(pheromone+ (node1*dim + node2)) = ((1 - alpha)*phi)+(alpha*delta);
  // *(pheromone+ (node2*dim + node1))= *(pheromone+ (node1*dim + node2));
}

int *two_opt(int *sol) {
  int best_gain = 1;
  int gain = 0;
  int i, j;
  int best_i, best_j = -1;
  while (best_gain != 0) {
    best_gain = 0;
    for (i = 0; i < dim - 1; i++) {
      for (j = i + 2; j < dim; j++) {
        gain = compute_gain(distances, sol, i, j);
        if (gain < best_gain) {
          best_gain = gain;
          best_i = i;
          best_j = j;
        }
      }
    }
    if (best_gain < 0) {
      two_opt_swap(sol, best_i, best_j);
      // exchange(sol, best_i, best_j);
      // check = check_valid(sol);
      // if (check == 1) {
      //     printf("%s\n", "invalid solution");
      //     exit(EXIT_FAILURE);
      // }
    }
  }
  return sol;
}

// wiki style
// int * two_opt(int * sol, int  * distances ){
//     int best_gain = 1;
//     int gain = 0;
//     int i, j;
//     int best_i, best_j = -1;
//     // int best_dist = calculate_solution_cost(sol, dim, distances);
//     // printf("%d\n", best_dist);
//     // int new_dist = 0;
//     int improvement = 1;
//     int * s = sol;
//     // for (size_t i = 0; i < dim; i++) {
//     //     printf("%d\t", s[i]);
//     // }
//     // printf("%s\n", "");
//
//     // printf("%d\n", calculate_solution_cost(s, dim, distances));
//
//     do {
//         for ( i = 0; i < dim -1; i++) {
//             for (j = i+1; j < dim; j++) {
//                 gain = compute_gain(distances, s, i, j);
//                 // two_opt_swap(s, i, k);
//                 // new_dist = calculate_solution_cost(s, dim, distances);
//                 // printf("%d\n", new_dist);
//                 if (gain < best_gain) {
//                     // printf("%d\t%d\n", new_dist, best_dist);
//                     best_gain  = gain;
//                     best_i = i;
//                     best_j = j;
//                 }
//                 else{
//                     improvement = 0;
//                 }
//             }
//         }
//         if (best_gain < 0) {
//             two_opt_swap(s, best_i, best_j);
//         }
//     } while(improvement == 1);
//     // for (size_t i = 0; i < dim; i++) {
//     //     printf("%d\t", s[i]);
//     // }
//     // printf("%s\n", "");
//     // printf("%d\n", 5/2);
//     // printf("%d\n", calculate_solution_cost(s, dim, distances));
//
// return s;
// }

void two_opt_swap(int *s, int i, int k) {
  int temp;
  int a = (i + 1) % dim;
  int b = k;
  while (a < b) {
    temp = s[a];
    s[a] = s[b];
    s[b] = temp;
    a = (a + 1) % dim;
    b = (b - 1) % dim;
  }
}
//
// int * hill_climbing(int * sol, int  * distances ){
//     int cost = 1;
//     int delta = 0;
//     int i, j;
//     int best_i, best_j = -1;
//     int stuck;
//     do {
//         stuck = 1;
//         for (i = 0; i < dim -1 ; i++) {
//             for (j = i+1; j < dim; j++) {
//                 delta = compute_gain(distances, sol ,i,j);
//                 if (delta < cost) {
//                     cost = delta;
//                     best_i = i;
//                     best_j = j;
//                     stuck = 0;
//                 }
//             }
//         }
//         if (cost < 0) {
//             two_opt_swap(sol, best_i, best_j);
//             // exchange(sol, best_i, best_j);
//         }
//     } while(stuck == 0);
//     return sol;
// }

//
// int * two_opt_first(int * sol, int  * distances ){
//     int best_gain = 1;
//     int gain = 0;
//     int i, j;
//     int best_i, best_j = -1;
//     while (best_gain !=0) {
//         best_gain = 0;
//         for (i = 0; i < dim -1 ; i++) {
//             for (j = i+1; j < dim; j++) {
//                 gain = compute_gain(distances, sol ,i,j);
//                 // printf("%s\t%d\n", "best gain", best_gain);
//                 // printf("%s\t%d\n", " gain", gain);
//
//                 // printf("%d\n", gain);
//                 if (gain < best_gain) {
//                     // printf("%s\n", " hello");
//                     best_gain = gain;
//                     best_i = i;
//                     best_j = j;
//                     break;
//                 }
//             }
//             if (best_gain < 0) {
//                 break;
//             }
//         }
//         if (best_gain < 0) {
//             exchange(sol, best_i, best_j);
//         }
//
//     }
//     return sol;
// }
int compute_gain(int *distances, int *sol, int i, int j) {
  // optimise this later
  // int a = i;
  // int b = (i + 1)%dim;
  // int c = j;
  // int d = (j + 1)%dim;
  if (i == j || i == (j + 1) % dim || (i + 1) % dim == j ||
      (i + 1) % dim == (j + 1) % dim) {
    return INT_MAX;
  }
  int gain = *(distances + (sol[i] * dim + sol[j])) +
             *(distances + (sol[(i + 1) % dim] * dim + sol[(j + 1) % dim]));
  int old_gain = *(distances + (sol[j] * dim + sol[(j + 1) % dim])) +
                 *(distances + (sol[i] * dim + sol[(i + 1) % dim]));
  return (gain - old_gain);
}

void exchange(int *sol, int best_i, int best_j) {
  for (size_t k = 0; k <= (best_j - best_i) / 2; k++) {
    swap_places(sol, best_i + k, best_j - k + 1);
  }
}
