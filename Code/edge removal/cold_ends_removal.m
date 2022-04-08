function [known_network, unknown_network] = cold_ends_removal(network, percent_removed)
%COLD_ENDS_REMOVAL This function removes edges with cold ends. It removes
%   nonbridge edges from low degree vertices.
%
%   NOTE: This is a deterministic function
%
%   Joshua Pickard jpic@umich.edu
%   April 8, 2022

    %% Set parameters
    known_edges = round(percent_removed * sum(sum(network)));
    current_edges = known_edges;
    
    %% Define the known edges
    known_network = network;
    bridges = false(size(network));
    
    %% Edge removal
    while current_edges >= known_edges
        % Select node with lowest degree
        [~, lowest_degree_node] = min(sum(known_network));
        % Select the neighbring node with the lowest degree
        neighbors = find(known_network(highest_degree_node,:) == true);
        [~, lowest_degree_neighbor] = min(sum(known_network(neighbors,:), 2));
        neighbor = neighbors(lowest_degree_neighbor);
        if ~isbridge(known_network, lowest_degree_node, neighbor)
            % Remove edge
            known_network(lowest_degree_node, neighbor) = false;
            known_network(neighbor, lowest_degree_node) = false;
            current_edges = current_edges - 2;
        else
            % Remove the edge for now (it will be added in later)
            known_network(lowest_degree_node, neighbor) = false;
            known_network(neighbor, lowest_degree_node) = false;
            % Count the edge as a bridge
            bridges(lowest_degree_node, neighbor) = true;
            bridges(neighbor, lowest_degree_node) = true;
        end
    end
    known_network = known_network + bridges;
    
    %% Define the unknown edges
    unknown_network = (network & ~known_network);
end

