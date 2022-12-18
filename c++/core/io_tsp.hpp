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
  bool exist_opt_ = false;
  std::vector<int> optimal_tour_;
  std::vector<std::vector<int>> dist_matrix_;
  std::vector<std::vector<double>> pheromone_;
  std::string name_;
  int nPoints_;
  int best_sol_;
  std::vector<std::vector<double>> points_;
  std::string file_name_;

  std::vector<int> new_sol;

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

    points_ =
        std::vector<std::vector<double>>(2, std::vector<double>(nPoints_));

    // store data

    for (int i = 7; i < nPoints_; i++) {
      boost::split(t, data[i], boost::is_any_of(" "));
      int index = stoi(t[0]);
      double x = stoi(t[1]);
      double y = stoi(t[2]);

      std::cout << index << std::endl;
      std::cout << x << std::endl;
      std::cout << y << std::endl;
      points_[0][index - 1] = x;
      points_[1][index - 1] = y;
    }

    create_dist_matrix();
    pheromone_ = std::vector<std::vector<double>>(
        nPoints_, std::vector<double>(nPoints_));

    for (int m = 0; m < nPoints_; m++) {
      for (int k = 0; k < nPoints_; k++) {
        pheromone_[m][k] = 1;
        pheromone_[k][m] = 1;
      }
    }

    new_sol = std::vector<int>(nPoints_);
  }

  void print_info() {
    std::cout << "\n\n#############################\n";
    std::cout << "name: " << name_ << std::endl;
    std::cout << "nPoints: " << nPoints_ << std::endl;
    std::cout << "best_sol: " << best_sol_ << std::endl;
    std::cout << "exist optimal: " << exist_opt_ << std::endl;
  }

  void create_dist_matrix() {
    dist_matrix_ =
        std::vector<std::vector<int>>(nPoints_, std::vector<int>(nPoints_));

    for (int i = 0; i < nPoints_; i++) {
      for (int j = 0; j < nPoints_; j++) {
        if (j == i) {
          dist_matrix_[i][j] = 9999;
          dist_matrix_[j][i] = 9999;
        } else {
          double x_sq = pow(points_[0][i] - points_[0][j], 2);
          double y_sq = pow(points_[1][i] - points_[1][j], 2);
          dist_matrix_[i][j] = sqrt(x_sq + y_sq);
          dist_matrix_[j][i] = sqrt(x_sq + y_sq);
        }
      }
    }
  }

  int getnPoints() { return nPoints_; }
};

#endif