function [theta, likelihoods] = NaiveKronFit(A, v, debug)
%NAIVEKRONFIT This function is slow but does a brute force evaluation of
% the kron fit problem according to equation 5.5 in Jure Leskovec's thesis.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

if nargin == 1
    v = false;
    debug = false;
elseif nargin == 2
    debug = false;
end

lr = 1e-5;

n = size(A,1);
k = length(size(A));

theta = rand(2 * ones(1, k));

maxItrs = 50;
likelihoods = zeros(maxItrs, 1);
for itr=1:maxItrs
    % Evaluate likelihood and gradient
    [l, gradients] = evaluateGradient(A, theta, debug);
    % Update model parameters
    theta = theta + lr * gradients;

    lr = 0.95 * lr;

    % Outputs
    if mod(itr, 5) == 0
        disp(itr);
        figure; plot(real(likelihoods));
    end
    likelihoods(itr) = l;
    if v 
        disp(l);
        disp(theta);
        disp(gradients * lr)
    end
end

end

