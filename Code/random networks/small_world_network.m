function [network, cords] = small_world_network(V, k, beta)
%SMALL_WORLD_NETWORK: This function generates a small world network. A
%   small world network is a network where most nodes are not neighbors,
%   but nodes with common neighbors are more likely to become neighbors. In
%   these networks, the average distance between any 2 nodes it
%   proportional to the Log of the number of nodes in the network. 
%
%   The Watts-Strogatz or Watts-Beta algorithm is used to construct small
%   world graphs in this implementation. This algorithm has input:
%   - V: the number of vertices in the graph
%   - k: the mean degree of the graph
%   - beta: a parameter indicating the randomness of the graph. It is
%           bounded by 0 <= beta <= 1 where 0 is a lattice and 1 is a
%           completely random graph
%   A general condition on these parameters is:
%                          N >> k >> ln(N) >> 1
%   N >> k insures 
%   1. Constructing a regular ring lattice is the first step.
%   2. All nodes with probability beta change one of their vertices to a
%   vertex outside its standard cluster
%
% Joshua Pickard jpic@umich.edu
% April 7, 2022

% 1. Construct a ring lattice
network = zeros(V,V);
for i=1:k % k is half of the average degree
    edges = ones(V, 1);
    inner = edges(1:length(edges)-i);       % Edges close to the diagonal
    outer = edges(length(edges)-i+1:end);   % Edges far from the diagonal
    network = network + diag(inner, i) + diag(outer, V-i); % Insert edges on off diagonals
end
network = logical(network);
network = network + network';

% Get cords to plot circular lattice
G = graph(network);
figure; 
p = plot(G, 'layout', 'circle');
x = p.XData;
y = p.YData;
cords = [x; y];
close;

% 2. Rewire all edges with probability beta
for v=1:V
    % For node v, select half of its edges (the set that may be rewired)
    possible_reconnect = zeros(V,1);
    if v+k <= V
        possible_reconnect(v+1:v+k) = true;
    else
        % Back connectiuons
        possible_reconnect(v+1:end) = true;
        back_connections = sum(possible_reconnect);
        front_connections = k-back_connections;
        possible_reconnect(1:front_connections) = true;
    end
    % Determine which of these need to be rewired
    rewire = possible_reconnect & rand(V,1) < beta;
    for v2=1:V
        if ~rewire(v2)
            continue
        end
        % Remove old edge
        network(v,v2) = false;
        network(v2,v) = false;
        % Select the new v2
        possible_neighbors = ~network(v,:);
        possible_neighbors(v) = false; % Avoid self loops
        new_neighbor = randsample(find(possible_neighbors == 1), 1);
        % Add new edge
        network(v,new_neighbor) = true;
        network(new_neighbor,v) = true;
    end
end
end

