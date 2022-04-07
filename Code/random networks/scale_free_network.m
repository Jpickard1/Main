function [network] = scale_free_network(V)
%scale_free_network: This function generates a scale free network using the
%   BA algorithm. A scale free network is a network where the degree
%   distribution of the nodes approximately follows the power law. This means
%   that the fraction of nodes with degree of at least k is proportional to
%   k^-g, where g is some constant.
%
%   The BA algorithm (after Albert-László Barabási and Réka Albert) uses
%   preferential attachment to construct scale free networks. The algorithm
%   can be broken down into the following stages:
%       1. Initialize a network with V0 nodes and E0 edges.
%       2. Loop until network is large enough
%           2.a Add a new node to the network
%           2.b Add edges between the new node and preexisting nodes where
%           the probability of an edge existing between the new node and
%           a preexisting node i is calculated as
%               P = degree of i / degree of all preexisting nodes
% Joshua Pickard jpic@umich.edu
% April 6, 2022

% 1. Intialize the network
network = false(V,V);
network(2,1) = true;
network(1,2) = true;

% 2. Loop until network is large enough
for v=3:V
    total_edges = sum(sum(network));
    for i=1:v
        if rand < sum(network(i, :)) / total_edges
            network(i,v) = true;
            network(v,i) = true;
        end
    end
end

end

