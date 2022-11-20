#ifndef UTILS_HPP
#define UTILS_HPP
#include <vector>

int compute_length(std::vector<int> solution, std::vector<std::vector<int> > dist_matrix){
    int total_length = 0;
    int starting_node = solution[0];
    int from_node = starting_node;
    for (int node : solution){
        total_length += dist_matrix[from_node][node];
        from_node = node;
    }

    total_length += dist_matrix[starting_node][from_node];
    return total_length;
}



#endif