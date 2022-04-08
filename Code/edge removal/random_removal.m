function [known_network, unknown_network] = random_removal(network, percent_removed)
    %% Set parameters
    known_edges = round(percent_removed * sum(sum(network)));
    current_edges = known_edges;

    %% Define the known edges
    known_network = network;
    [m, ~] = size(network);
    while current_edges >= known_edges
        for i=1:m
            for j=i:m
                if current_edges >= known_edges
                    break
                end
                if network(i,j) && rand() > prob
                    known_network(i,j) = 0;
                    known_network(j,i) = 0;
                    current_edges = current_edges - 2;
                end
            end
        end
    end

    %% Define the unknown edges
    unknown_network = (network & ~known_network);
end
