function [network] = erdos_renyi_gilbert_network(V, p)
%erdos_renyi_gilbert_network: This function generates an Erdos-Renyi-Gilbert
%   random network. An Erdos-Renyi-Gilbert graph is a graph with V vertices
%   where every edge occurs with probability p.
%
% Joshua Pickard jpic@umich.edu
% April 6, 2022

edges = (rand(nchoosek(V, 2), 1) < p);
network = logical(squareform(edges));

end
