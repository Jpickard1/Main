%% Evaluate Gradient
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [likelihood, gradient]=evaluateGradient(A, theta, debug)

% TODO: rest
itrs = 100;

n0 = size(theta,1);
n = size(A,1);
kronExp = log(n) / log(n0);

% precalculate likeliehood of empty graph
le = getEmptyGraphLL(n, theta);

% precalculate the gradient of empty graph
ge = getEmptyGraphGrad(n, theta);

% precompute edge set as a matrix
[e1, e2] = find(A ~= 0);
E = [e1 e2];

likelihoods = zeros(itrs, 1);
% gradient    = zeros(size(theta));
gradients   = cell(itrs, 1);
for t=1:itrs
    % Generate sample permutation
    if t == 1
        p = samplePermutation(A, theta);
    else
        p = samplePermutation(A, theta, p);
    end
    % if debug; disp(p); end

    % Calculate log likelihood
    likelihood = le;
    for e=1:size(E,1)
        eLL = edgeLL(n, theta, p(E(e,1)), p(E(e,2)));
        likelihood = likelihood - log(1 - exp(eLL)) + eLL;
    end
    likelihoods(t) = likelihood;

    % Calculate gradient
    gradUpdate = ge;
    for e=1:size(E,1)
        eGrad = edgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
        gradUpdate = gradUpdate + eGrad;
        % gradUpdate = gradUpdate - log(1 - eGrad) + log(eGrad);
    end
    gradients{t} = real(gradUpdate);
end

likelihood = mean(likelihoods);
gradient = zeros(size(theta));
for i=1:itrs
    gradient = gradient + gradients{i};
end
gradient = gradient ./ itrs;

end
