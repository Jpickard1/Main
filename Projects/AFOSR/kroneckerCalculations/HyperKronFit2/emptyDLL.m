function [DLL] = emptyDLL(n, theta)
%EMPTYLL Computes the derivative of the log likelihood of an empty 
% Kronecker graph with respect to theta

n0 = size(theta,1);

kronExp = ceil(log(n) / log(n0));

s  = sum(theta, 'all');
sq = sum(theta .^ 2, 'all');

DLL = -kronExp * (s ^ (kronExp - 1)) - kronExp * (sq ^ (kronExp - 1)) * theta;

end