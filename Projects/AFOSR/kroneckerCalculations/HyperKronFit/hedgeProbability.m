function [P] = hedgeProbability(n, theta, idxs)
%HEDGELLPROBABILITY Evaluates the probability of a hyperedge
%
% PARAMETERS
%   n is number of vertices in a hypergraph
%   theta is the kronecker initiator parameters
%   idxs are the vertices of the hyperedge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

theta = theta / sum(sum(theta));

n0 = size(theta,1);
kronExp = log(n) / log(n0);

P = 1;
for i=1:kronExp
    idx = mod(floow((idxs - 1) / n0^(i - 1)), n0) + 1;
    idx = cell2mat(idx);
    P = P * theta(idx{:});
end


end