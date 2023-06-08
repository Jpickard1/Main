function [theta, likelihoods] = NaiveKronFit(A, v, debug, n0)
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

% theta = rand(2 * ones(1, k));
theta = [0.9 0.6;
         0.6 0.1];

if nargin == 4
    theta = rand(n0, n0);
end
n0 = size(theta,1);

maxItrs = 25;
likelihoods = zeros(maxItrs, 1);
for itr=1:maxItrs
    if v
        fprintf("Itr: %d ]\n", itr);
    end
    
    % Evaluate likelihood and gradient
    % [l, gradients] = evaluateGradient(A, theta, debug);
    [l, gradients] = sampleGradient(A, theta, debug);
    % Update model parameters
    thetaOld = theta;
    theta = theta + lr * gradients;

    lr = 0.95 * lr;
    for i=1:n0
        for j=1:n0
            if theta(i,j) > 0.9999; theta(i,j) = 0.9999; end
            if theta(i,j) < 0.0; theta(i,j) = 0.0; end
        end
    end

    if v
        fprintf("CurrentLL: %d\n", l);
        fprintf("Gradient Updates:\n");
        for i=1:numel(theta)
            fprintf("    %d]  %f = %f + %f\n", i, theta(i), thetaOld(i), lr * gradients(i));
        end
        fprintf(' \n');
    end

    % Outputs
    if mod(itr, 5) == 0
        disp(itr);
        figure; plot(real(likelihoods));
        pause(1);
    end
    likelihoods(itr) = l;
end

end

