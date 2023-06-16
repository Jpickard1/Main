function LL = getEmptyHypergraphLLapx(n, theta)
%GETEMPTYGRAPHGRAD Returns the log likelihood of a hypergraph with 0 edges.
%
%   The log likelihood (LL) of an edge not existing is log(1-x), which has
%   a second order Taylor expansion:
%
%       log(1-x) ~ -x -0.5^x
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

theta = theta / sum(theta, 'all');

n0 = size(theta,1);
kronExp = log(n) / log(n0);

s  = sum(theta,'all');
sq = sum(theta .^ 2,'all');

LL = -1 * (s^kronExp) - 0.5 * (sq ^ kronExp);

end
