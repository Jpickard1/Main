function [network] = erdos_renyi_network(V, E)
%erdos_renyi_network: This function generates an erdos-renyi random network
%   An Erdos-Renyi graph is a graph with V vertices and E edges where the
%   graph is chosen uniformly at random from the set of all graphs with V
%   vertices and E edges. The vertices are ordered in the sense that graphs
%   that can be obtained from one another via graph permutations are
%   distinct.
%
% Joshua Pickard jpic@umich.edu
% April 6, 2022

num_edges_total = nchoosek(V, 2);
false_edges = zeros(num_edges_total - E, 1);
true_edges = ones(E, 1);
ordered_edges = [false_edges; true_edges];
shuffled_edges = ordered_edges(randperm(length(ordered_edges)));
network = logical(squareform(shuffled_edges));
end

