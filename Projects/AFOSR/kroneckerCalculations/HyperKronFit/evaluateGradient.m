%% Evaluate Gradient
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [likelihood, gradient]=evaluateGradient(A, theta, debug)

if nargin == 2
    debug = false;
end

% TODO: rest
itrs = 10;

n0 = size(theta,1);
n = size(A,1);
k = log(n) / log(n0);
n1 = size(theta,1);

% precalculate likeliehood of empty graph
le = -1 * (sum(sum(theta)) ^ k);
s1 = 0;
for i=1:size(theta,1)
    for j=1:size(theta,2)
        s1 = s1 + theta(i,j)^2;
    end
end
le = le + 0.5 * s1^k;

% precalculate the gradient of empty graph
ge = -1 * ones(n1, n1) - theta;

% precompute edge set as a matrix
[e1, e2] = find(A ~= 0);
E = [e1 e2];

likelihoods = zeros(itrs, 1);
gradient    = zeros(size(theta));
for t=1:itrs
    % Generate sample permutation
    p = samplePermutation(A, theta);
    if debug; disp(p); end

    % Calculate log likelihood
    likelihood = le;
    for e=1:size(E,1)
        eProb = log(edgeProbability(n, theta, p(E(e,1)), p(E(e,2))));
        likelihood = likelihood - log(1 - eProb) + log(eProb);
    end
    likelihoods(t) = likelihood;

    % Calculate gradient
    gradUpdate = ge;
    for e=1:size(E,1)
        eGrad = edgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
        gradUpdate = gradUpdate - log(1 - eGrad) + log(eGrad);
    end
    gradient = gradient + real(gradUpdate);
end

likelihood = mean(likelihoods);
gradient = gradient ./ itrs;

end
