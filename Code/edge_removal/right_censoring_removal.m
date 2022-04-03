function [known_network, unknown_network] = right_censoring_removal(network)
    %% Set parameters
    max_degree = max(sum(network));
    p_degree = 0.8;
    max_allowed_degree = round(max_degree * p_degree);
    
    %% Define the known edges
    known_network = network;
    
    %% Edge removal
    while max(sum(known_network)) > max_allowed_degree
        % Select node with highest degree
        [~, highest_degree_node] = max(sum(known_network));
        % Select the neighbring node with the highest degree
        neighbors = find(known_network(highest_degree_node,:) == true);
        % neighbors = [1 3 5 7 9]
        [~, highest_degree_neighbor] = max(sum(known_network(neighbors,:), 2));
        neighbor = neighbors(highest_degree_neighbor);
        % Remove edge
        known_network(highest_degree_node, neighbor) = false;
        known_network(neighbor, highest_degree_node) = false;
    end
    
    %% Define the unknown edges
    unknown_network = (network & ~known_network);
end 
