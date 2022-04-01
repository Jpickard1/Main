function [known_network, unknown_network] = random_removal(network)
    %% Set parameters
    prob = 0.5; % Probability of an edge being known

    %% Define the known edges
    known_network = network;
    [m, ~] = size(network);
    for i=1:m
        for j=i:m
            if network(i,j) && rand() > prob
                known_network(i,j) = 0;
                known_network(j,i) = 0;
            end
        end
    end

    %% Define the unknown edges
    unknown_network = (network & ~known_network);
end
