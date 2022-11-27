#ifndef IO_TSP_HPP
#define IO_TSP_HPP
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>

class ProblemInstance {
public:
  std::vector<std::string> read_file(std::string file_name) {
    std::string line;
    std::vector<std::string> data;
    std::ifstream mFile(file_name);
    while (std::getline(mFile, line)) {
      data.push_back(std::move(line));
    }
    mFile.close();

    return data;
  }

  ProblemInstance(const std::string name_tsp) {
    file_name_ = name_tsp;
    std::vector<std::string> data = read_file(file_name_);
  }

private:
  bool exist_opt_ = false;
  std::vector<int> optimal_tour_;
  std::vector<std::vector<int>> dist_matrix_;
  std::string name_;
  std::vector<int> nPoints_;
  std::vector<float> best_sol_;
  std::string file_name_;
};

#endif