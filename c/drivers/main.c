// #include "acs.h"
#include "io_tsp.h"
// #include "utils .h"

int main(int argc, char const *argv[]) {

  FILE *fp;

  // Read file
  if (argc < 2) {
    printf("%s\n", "please provide the filepath for the tsp.");
    exit(EXIT_FAILURE);
  } else if (argc == 2) {
    fp = fopen(argv[1], "r");
    if (!fp) {
      printf("%s\n", "invalid filepath");
      exit(EXIT_FAILURE);
    }
  } else {
    return 1;
  }

  IO *io;
  read_words(fp, io);
  printf("%i", io->dim);
  //   int dim = io.get_dim();
  //   double nodes[2][dim];
  //   int distances[dim][dim];
  //   double pheromone[dim][dim];
  //   int sol[dim];
  //   int new_sol[dim];

  //   io.read_coords(fp, nodes);
  //   calculate_distance(nodes, distances);

  //   // Setup pheromone matrix.
  //   for (size_t m = 0; m < dim; m++) {
  //     for (size_t k = 0; k < dim; k++) {
  //       pheromone[m][k] = 1;
  //       pheromone[k][m] = 1;
  //     }
  //   }

  //   // Populate solution array with 1's.
  //   for (size_t i = 0; i < dim; i++) {
  //     sol[i] = i;
  //   }

  //   // srand(time(NULL));

  //   // Something
  //   for (size_t j = dim - 1; j > 0; j--) {
  //     int r = rand() % dim;
  //     swap_places(sol, j, r);
  //   }

  //   printf("%s\t%d\n", "Best: ", best);

  //   /*
  //    * Make ants
  //    */

  //   int swarm_count = 3;
  //   ANT *nest[swarm_count];
  //   ANT *a1;
  //   int test = 1;
  //   for (int i = 0; i < swarm_count; i++) {
  //     // initialise ant
  //     nest[i] = (ANT *)malloc(sizeof(ANT));
  //     a1 = nest[i];
  //     a1->trail = (int *)malloc(dim * sizeof(int));
  //     a1->visited = (int *)malloc(dim * sizeof(int));
  //   }
  //   while (test != 0) {
  //     acs(nest, swarm_count, distances, pheromone, new_sol);
  //   }

  //   int sol_cost = calculate_solution_cost(new_sol, sizeof(new_sol) /
  //   sizeof(int),
  //                                          distances);
  //   printf("%s\t%d\n", "ACS: ", sol_cost);

  //   for (int i = 0; i < swarm_count; i++) {
  //     free(nest[i]->trail);
  //     free(nest[i]->visited);
  //     free(nest[i]);
  //   }
  //   // nearest_neighbor(sol, distances,new_sol);
  //   // sol_cost = calculate_solution_cost(new_sol,
  //   sizeof(new_sol)/sizeof(int),
  //   // distances); printf("%s\t%d\n","NN: " ,sol_cost);
  //   printf("%s\t%f\n", "Err: ", ((double)(sol_cost - best) / (double)best) *
  //   100); int valid; valid = check_valid(new_sol); if (valid == 0) {
  //     FILE *fw = fopen("near.txt", "w");
  //     if (fw) {
  //       for (size_t k = 0; k < dim; k++) {
  //         // printf("%d\t", new_sol[k]);
  //         fprintf(fw, "%d\n", new_sol[k] + 1);
  //       }
  //       fclose(fw);
  //     }
  //   }
  //   fclose(fp);
  // }
  // return 0;
}