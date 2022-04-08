function [network] = small_world_network(V, k, beta)
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
%   1. Constructing a regular lattice is the first step. The V edges are
%   connected to the k vertices before and the k vertices after it (where
%   the vertices are arranged in a circle)
%   2. All nodes with probability beta change one of their vertices to a
%   vertex outside its standard cluster
%
% Joshua Pickard jpic@umich.edu
% April 7, 2022

V = 4;
k = 1;
s = repelem((1:V)',1,k);
t = s + repmat(1:k,V,1);
t = mod(t-1,V)+1;
h = graph(s,t);
plot(h)

clear
V = 5;
A = diag(ones(V-1,1),1)+diag(ones(V-2,1),2)+diag(ones(2,1),V-2)+diag(ones(1),V-1);

end

