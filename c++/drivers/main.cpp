#include "acs.hpp"
#include "io_tsp.hpp"
#include "utils.hpp"
#include <filesystem>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;

void run(const bool show_plots = false, const bool verbose = false) {
  fs::path parent_folder = fs::current_path().parent_path().parent_path();
  fs::path folder =
      parent_folder.string() + "/AI_cup_2022/" + "AI_cup_2022_problems/";
  std::vector<std::string> file_list;

  if (!std::filesystem::is_directory(folder)) {
    throw std::runtime_error(folder.string() + " is not a folder");
  }

  for (const auto &entry : std::filesystem::directory_iterator(folder)) {
    fs::path full_name = entry.path();
    full_name = fs::absolute(full_name);
    full_name = full_name.string();

    if (entry.is_regular_file()) {
      const auto base_name = entry.path().filename().string();
      file_list.push_back(full_name);
    }
  }

  // Iterate over the problems

  ProblemInstance problem(file_list[0]);
  problem.print_info();

  // for (auto problem_path : file_list) {
  //   // Create ProblemInstance.
  //   ProblemInstance problem(problem_path);
  //   // Choose solver.
  //   problem.print_info();

  //   // Chose improve.
  // }

  int swarm_count = 3;
  std::vector<ANT> nest(swarm_count);
  ANT a1;
  int test = 1;
  // for (int i = 0; i < swarm_count; i++) {
  //   // initialise ant
  //   // nest[i] = (ANT *)malloc(sizeof(ANT));
  //   // a1 = nest[i];
  //   a1.trail = (int *)malloc(problem.nPoints_ * sizeof(int));
  //   a1.visited = (int *)malloc(problem.nPoints_ * sizeof(int));
  // }
  while (test != 0) {
    acs(nest, swarm_count, problem.dist_matrix_, problem.pheromone_,
        problem.new_sol, problem.nPoints_);
  }

  // int sol_cost = calculate_solution_cost(new_sol, sizeof(new_sol) /
  // sizeof(int),
  //                                        distances);
  // printf("%s\t%d\n", "ACS: ", sol_cost);

  // for (int i = 0; i < swarm_count; i++) {
  //   free(nest[i]->trail);
  //   free(nest[i]->visited);
  //   free(nest[i]);
  // }
  // // nearest_neighbor(sol, distances,new_sol);
  // // sol_cost = calculate_solution_cost(new_sol, sizeof(new_sol)/sizeof(int),
  // // distances); printf("%s\t%d\n","NN: " ,sol_cost);
  // printf("%s\t%f\n", "Err: ", ((double)(sol_cost - best) / (double)best) *
  // 100); int valid; valid = check_valid(new_sol); if (valid == 0) {
  //   FILE *fw = fopen("near.txt", "w");
  //   if (fw) {
  //     for (size_t k = 0; k < dim; k++) {
  //       // printf("%d\t", new_sol[k]);
  //       fprintf(fw, "%d\n", new_sol[k] + 1);
  //     }
  //     fclose(fw);
  //   }
  // }
  // fclose(fp);
}

int main(int argv, const char *argc[]) { run(); }