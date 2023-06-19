function [nDLL] = noHEdgeDLL(n, theta, idxs)
%HEDGE evaluates the log likelihood of an edge not existing in a hypergraph
%
%   (d / d theta_i,j) log(1 - LL) =
%   (d / d theta_i,j) log(1 - theta(i,j)^c(i,j)*theta(k,l)^c(k,l) ... ) = 
%       (c(i,j) * LL / theta(i,j)) / (1 - LL)
%       
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

LL = hedgeLL(n, theta, idxs);
counts = getCounts(n, theta, idxs);

LL = LL * ones(size(theta));
nDLL = (counts * LL ./ theta) ./ (1 - LL);

end