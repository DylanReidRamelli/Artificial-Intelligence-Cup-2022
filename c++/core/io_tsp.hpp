#ifndef IO_TSP_HPP
#define IO_TSP_HPP
#include "utils.hpp"
#include <boost/algorithm/string.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
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

    // store data set information
    std::vector<std::string> t;
    boost::split(t, data[0], boost::is_any_of(": "));
    name_ = t[t.size() - 1];
    boost::split(t, data[3], boost::is_any_of(": "));
    nPoints_ = std::stoi(t[t.size() - 1]);
    boost::split(t, data[5], boost::is_any_of(": "));
    best_sol_ = std::stoi(t[t.size() - 1]);

    create_dist_matrix();
  }

  void print_info() {
    std::cout << "\n\n#############################\n";
    std::cout << "name: " << name_ << std::endl;
    std::cout << "nPoints: " << nPoints_ << std::endl;
    std::cout << "best_sol: " << best_sol_ << std::endl;
    std::cout << "exist optimal: " << exist_opt_ << std::endl;
  }

  void create_dist_matrix() {
    dist_matrix_.resize(nPoints_ * nPoints_);
    for (int i = 0; i < nPoints_; i++) {
      for (int j = i; j < nPoints_; j++) {
        std::vector<int> rowi = points_[i];
        std::vector<int> rowj = points_[j];
        auto firsti = rowi.begin() + 1;
        auto lasti = rowi.begin() + 3;
        auto firstj = rowj.begin() + 1;
        auto lastj = rowj.begin() + 3;
        std::vector<int> subi(firsti, lasti);
        std::vector<int> subj(firstj, lastj);

        dist_matrix_[j][i] = distance_euc(subi, subj);
      }
    }
  }

private:
  bool exist_opt_ = false;
  std::vector<int> optimal_tour_;
  std::vector<std::vector<int>> dist_matrix_;
  std::string name_;
  int nPoints_;
  int best_sol_;
  std::vector<std::vector<int>> points_;
  std::string file_name_;
};

#endif