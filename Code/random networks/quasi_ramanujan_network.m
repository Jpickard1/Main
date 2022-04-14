function [network, cords] = quasi_ramanujan_network(V, k, beta)
%QUASI_RAMANUJAN_NETWORK This function generates a random, quasi ramanujan
%   network based on the paper "Algebraic Connectivity Rtio of Ramanujan 
%   Graphs" by Reza Olifati-Saber. The paper proposes a method for
%   generating these networks by first creating a Small World network using
%   the Watts-Strogatz method and then balancing the network to be d
%   regular with Algorithm 1 in the paper.
beta = 1;
[network, cords] = small_world_network(V, k, beta);

%sum(network);
%i = 0;
while max(sum(network)) ~= min(sum(network))
    %disp("iter: " + string(i))
    %i = i + 1;
    %disp(max(sum(network)))
    %disp(min(sum(network)))
    % Select node with max degree
    rich_node = randsample(find(sum(network) == max(sum(network))), 1);
    % Select node with min degree
    poor_node = randsample(find(sum(network) == min(sum(network))), 1);
    % Select random neighbor of rich node
    rich_neighbor = randsample(find(network(rich_node,:) == 1), 1);
    % Rewire rich_node <--> rich_neighbor to rich_neighbor <--> poor_node
    network(rich_node, rich_neighbor) = false;
    network(rich_neighbor, rich_node) = false;
    network(poor_node, rich_neighbor) = true;
    network(rich_neighbor, poor_node) = true;
end
%disp("iter: " + string(i))
%i = i + 1;
%disp(max(sum(network)))
%disp(min(sum(network)))

end

