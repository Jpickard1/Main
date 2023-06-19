function [LL] = hedgeLL(n, theta, idxs)
%HEDGE evaluates the log likelihood of an edge existing in a hypergraph
%
% ARGUMENTS
%   n is size of Kronecker hypergraph
%   theta is Kronecker inidiator
%   idxs is the node in the Kronecker graph
%
% TODO
%   this is only configured for graphs. see noHEdgeLL.m
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

n0 = size(theta,1);
counts = getCounts(n, theta, idxs);
lTheta = log(theta);
LL = sum(counts .* lTheta, 'all');

end