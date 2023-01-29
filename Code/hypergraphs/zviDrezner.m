function multirelation = zviDrezner(data)
% ZVI DREZNER
%   This function computes the multiway relationship between random
%   variables of data according to the formula given by Zvi Drezner in the
%   paper: "Multirelation — a correlation among more than two variables"
%
%   data: m x n matrix with m measurements on n random variables
%
%   Ref: https://www.sciencedirect.com/science/article/pii/0167947393E00467
%
% Auth: Joshua Pickard (jpic@umich.edu)
% Date: September 14, 2022

C = corrcoef(data);
C(isnan(C)) = 0;
multirelation = 1 - min(eig(C));

end
