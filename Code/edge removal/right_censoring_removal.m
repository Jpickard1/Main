function [known_network, unknown_network] = right_censoring_removal(network, percent_removed)
%   NOTE: This is a deterministic function
    %% Set parameters
    known_edges = round(percent_removed * sum(sum(network)));
    current_edges = known_edges;
    
    %% Define the known edges
    known_network = network;
    
    %% Edge removal
    while current_edges >= known_edges
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
        current_edges = current_edges - 2;
    end
    
    %% Define the unknown edges
    unknown_network = (network & ~known_network);
end 
