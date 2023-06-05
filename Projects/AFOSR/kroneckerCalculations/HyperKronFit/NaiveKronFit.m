function [theta] = NaiveKronFit(A)
%NAIVEKRONFIT This function is slow but does a brute force evaluation of
% the kron fit problem according to equation 5.5 in Jure Leskovec's thesis.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

lr = 0.01;

n = size(A,1);
k = length(size(A));

theta = rand(2 * ones(1, k));

maxItrs = 10;
for itr=1:maxItrs
    disp(itr);
    [~, gradients] = evaluateGradient(A, theta);
    theta = theta + lr * sum(gradients);
end

end

