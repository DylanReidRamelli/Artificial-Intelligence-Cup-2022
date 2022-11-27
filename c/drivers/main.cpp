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

  // for (auto problem_path : file_list) {
  // Create ProblemInstance.
  ProblemInstance problem(file_list[0]);
  // std::cout << problem_path << std::endl;
  // Choose solver.

  // Chose improve.
  // }
}

int main(int argv, const char *argc[]) {

  std::vector<int> solution[4];
  std::vector<std::vector<int>> dist_matrix[4][4];

  run();

  // compute_length();
}