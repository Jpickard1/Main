function [LL] = emptyLL(n, theta)
%EMPTYLL Computes the log likelihood of an empty Kronecker graph

n0 = size(theta,1);
kronExp = log(n) / log(n0);

s  = sum(sum(theta));
sq = sum(sum(theta .^ 2));

LL = (-1 * (s^kronExp)) - (0.5 * (sq ^ kronExp));

end