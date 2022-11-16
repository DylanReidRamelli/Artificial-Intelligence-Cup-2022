def solution_length(solution, adjacency_matrix):
    solution_length = 0
    start_node = solution[0]
    for node in solution[1:]:
        solution_length += adjacency_matrix[start_node,node]
        start_node = node
    # Add edge from start node to last node.
    solution_length += adjacency_matrix[solution[0],start_node]    
    return solution_length