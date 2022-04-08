function [network] = scale_free_network(V, m0)
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
%
% Joshua Pickard jpic@umich.edu
% April 6, 2022

% 1. Intialize the network randomly
network = false(V,V);
m0_network = erdos_renyi_network(m0, round(2*m0));
network(1:m0,1:m0) = m0_network;

% 2. Loop until network is large enough
for v=m0+1:V
    total_edges = sum(sum(network));
    while sum(network(v,:)) == 0
        for i=1:v-1
            if rand < sum(network(i, :)) / total_edges
                network(i,v) = true;
                network(v,i) = true;
            end
            if sum(network(v,:)) >= m0
                break
            end
        end
    end
end

end

