function [known_network, unknown_network] = snowball_removal(network, percent_removed)
%   NOTE: This is a deterministic function
    %% Set parameters
    known_edges = round(percent_removed * sum(sum(network)));

    %% Define the known edges
    known_network = zeros(size(network));
    current_edges = 0;
    
    %% BFS to select the known edges
    %node = round(rand() * length(network));
    queue = [round(rand() * length(network))];
    known_nodes = [];
    while current_edges < known_edges
        % Set node and pop queue
        node = queue(1);
        queue = queue(2:end);
        % Skip if the node has already been discovered
        if any(known_nodes == node)
            continue;
        end
        % Find edges from the current node
        connections = network(node, :);
        new_nodes = find(connections == true);
        % Add known connections to the network
        for n=new_nodes
            known_network(n, node) = true;
            known_network(node, n) = true;
            current_edges = current_edges + 2;
        end
        % Update the queu and known nodes
        queue = [queue new_nodes];
        known_nodes = [known_nodes node];
    end 
    %% Define the unknown edges
    unknown_network = (network & ~known_network);
end
