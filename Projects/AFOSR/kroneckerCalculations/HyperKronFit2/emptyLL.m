function [LL] = emptyLL(n, theta)
%EMPTYLL Approximates the log likelihood of an empty Kronecker graph

n0 = size(theta,1);
kronExp = log(n) / log(n0);

s  = sum(theta, 'all');
sq = sum(theta .^ 2, 'all');

LL = (-1 * (s^kronExp)) - (0.5 * (sq ^ kronExp));

end