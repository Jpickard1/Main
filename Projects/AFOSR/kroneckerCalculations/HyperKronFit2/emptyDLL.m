function [DLL] = emptyDLL(n, theta)
%EMPTYLL Computes the derivative of the log likelihood of an empty 
% Kronecker graph with respect to theta

n0 = size(theta,1);

kronExp = log(n) / log(n0);

s  = sum(sum(theta));
sq = sum(sum(theta .^ 2));

DLL = -kronExp * (s ^ (kronExp - 1)) - kronExp * (sq ^ (kronExp - 1)) * theta;

end