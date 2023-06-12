function [LL] = hedgeLLapx(n, theta, idxs)
%HEDGELLAPX Approximates the log likelihood of a hyperedge
%
% PARAMETERS
%   n is number of vertices in a hypergraph
%   theta is the kronecker initiator parameters
%   idxs are the vertices of the hyperedge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

theta = theta / sum(theta, 'all');

n0 = size(theta,1);
kronExp = log(n) / log(n0);

LL = 0;
for i=1:kronExp
    idx = mod(floor((idxs - 1) / n0^(i - 1)), n0) + 1;
    idx = num2cell(idx);
    LL = LL - theta(idx{:}) - 0.5 * theta(idx{:})^2;
end


end