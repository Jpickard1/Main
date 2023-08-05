function [s] = numSolns(n, k)
%NUMSOLNS gives the number of solutions to tensor characteristic polynomial
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2023

s = (((k-1)^n) - 1) / (k-2);

end

