function multirelation = benjaminTaylor(data)
% BENJAMIN TAYLOR
%   This function computes the multiway relationship between random
%   variables of data according to the formula given by Benjamin Taylor in
%   the paper: "A Multi-Way Correlation Coefficient"
%
%   data: m x n matrix with m measurements on n random variables
%
%   Ref: https://arxiv.org/abs/2003.02561
%
% Auth: Joshua Pickard (jpic@umich.edu)
% Date: September 14, 2022

C = corrcoef(data);
multirelation = (1/sqrt(length(C))) * std(eig(C));

end
