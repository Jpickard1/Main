function multirelation = wangZheng(data)
% Wang Zheng
%   This function computes the multiway relationship between random
%   variables of data according to the formula given by Jianji Wang and 
%   Nanning Zheng in their paper: "Measures of Correlation for Multiple 
%   Variables."
%
%   data: m x n matrix with m measurements on n random variables
%
%   Ref: https://arxiv.org/abs/1401.4827
%
% Auth: Joshua Pickard (jpic@umich.edu)
% Date: September 14, 2022

C = corrcoef(data);
multirelation = (1 - det(C)) ^ 0.5;

end
