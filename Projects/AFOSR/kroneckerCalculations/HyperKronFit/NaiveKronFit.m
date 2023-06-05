function [theta, likelihoods] = NaiveKronFit(A, debug)
%NAIVEKRONFIT This function is slow but does a brute force evaluation of
% the kron fit problem according to equation 5.5 in Jure Leskovec's thesis.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

if nargin == 1
    debug = false;
end

lr = 0.001;

n = size(A,1);
k = length(size(A));

theta = rand(2 * ones(1, k));

maxItrs = 1000;
likelihoods = zeros(maxItrs, 1);
for itr=1:maxItrs
    % Evaluate likelihood and gradient
    [l, gradients] = evaluateGradient(A, theta, debug);
    % Update model parameters
    theta = theta + lr * gradients;

    % Outputs
    if mod(itr, 10) == 0; disp(itr); end
    likelihoods(itr) = l;
end

end

