function [likelihood, gradient]=sampleGradientH(A, theta, debug)
% SAMPLEGRADIENTH The hypergraph implementation of sample gradient.
%
%   This function aproximates the gradient of each parameter in the Kronecker
%   initiator tensor.
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023
 
% TODO: rest
itrs = 5000;

n0 = size(theta,1);
n = size(A,1);
k = length(size(theta));
kronExp = log(n) / log(n0);

% precompute edge set as a matrix
idxs = cell(1, k);
linIdx = find(A > 0);
[idxs{:}] = ind2sub(size(A), linIdx);
E = cell2mat(idxs);

p = firstPermutation(A, theta);                                            % This needs to be written 
if debug; disp(p); end
likelihoods = zeros(itrs, 1);
gradients   = cell(itrs, 1);

% Calculations for first permutation
le = getEmptyHypergraphLL(n, theta);
ge = getEmptyHypergraphGrad(n, theta);

% Calculate log likelihood
likelihood = le;
for e=1:size(E,1)
    eLL = hedgeLL(n, theta, p(E(e,1)), p(E(e,2)));
    likelihood = likelihood - log(1 - exp(eLL)) + eLL;
end
% Calculate gradient
gradUpdate = ge;
for e=1:size(E,1)
    eGrad = hedgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
    gradUpdate = gradUpdate + eGrad;
end
likelihoods(1) = likelihood;
gradients{1} = gradUpdate;

accepted = 0;
for t=2:itrs
    % if debug; disp(p); end
    % Generate sample permutation
    [p, accept, u, v] = nextPermutation(A, theta, p);

    % Update 
    if accept
        % Update number of accepted permutations
        accepted = accepted + 1;

        % Update gradient
        % gradUpdate = updateGradForPerm(A, theta, gradients{t - 1}, u, v);
        % Update likelihood
        % likelihood = updateLLForPerm(A, theta, likelihoods{t - 1}, u, v);

            % This section will be improved to the above commented out code
            % to improve the time performance.

            % Calculate log likelihood
            likelihood = le;
            for e=1:size(E,1)
                eLL = edgeLL(n, theta, p(E(e,1)), p(E(e,2)));
                likelihood = likelihood - log(1 - exp(eLL)) + eLL;
            end
            likelihoods(t) = likelihood;
        
            % Calculate gradient
            %   I could speed this loop up by summing only over the sets of
            %   edges that have changed (i.e. edges with endpoints u and v)
            gradUpdate = ge;
            for e=1:size(E,1)
                eGrad = hedgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
                gradUpdate = gradUpdate + eGrad;
            end
            gradients{t} = real(gradUpdate);

        % Save updated calculations
        likelihoods(t) = likelihood;
        gradients{t} = real(gradUpdate);
    else
        likelihoods(t) = likelihoods(t-1);
        gradients{t} = gradients{t-1};
    end
end

fprintf("Accepted: %d / %d \n", accepted, itrs);

likelihood = mean(likelihoods);
gradient = zeros(size(theta));
for i=1:itrs
    gradient = gradient + gradients{i};
end
gradient = gradient ./ itrs;
if debug; disp(p); end
end
