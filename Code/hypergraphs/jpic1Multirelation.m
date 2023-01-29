function multirelation = jpic1Multicorrelation(data)
% Jpic1 Multicorrelation
%   This function gives the multicorrelation between random variables in a
%   data matrix loosely based on the Tucker Decomposition (very loosely).
%   Rather than determining the multiway relationship between a series of
%   variables 
%
%   data: m x n matrix with m measurements on n random variables
%
% Auth: Joshua Pickard (jpic@umich.edu)
% Date: September 14, 2022

p = prod(data, 2);
multirelation = sum(p);

end
