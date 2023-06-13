function LL = getEmptyGraphLLapx(n, theta)
%GETEMPTYGRAPHGRAD Returns the log likelihood of a graph with 0 edges.
%
%   The log likelihood (LL) of an edge not existing is log(1-x), which has
%   a second order Taylor expansion:
%
%       log(1-x) ~ -x -0.5^x
%
%   SNAP: kronecker.cpp 987-994 GetApxEmptyGraphLL
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 6, 2023

% theta = theta / sum(sum(theta));

n0 = size(theta,1);
kronExp = log(n) / log(n0);

s  = sum(sum(theta));
sq = sum(sum(theta .^ 2));

LL = (-1 * (s^kronExp)) - (0.5 * (sq ^ kronExp));

end
