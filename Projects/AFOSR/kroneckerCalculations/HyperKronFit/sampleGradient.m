function [likelihood, gradient]=sampleGradient(A, theta, debug)
% SAMPLEGRADIENT
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023
 
% TODO: rest
itrs = 1000;

n0 = size(theta,1);
n = size(A,1);
kronExp = log(n) / log(n0);

% precompute edge set as a matrix
[e1, e2] = find(A ~= 0);
E = [e1 e2];

p = firstPermutation(A, theta);
likelihoods = zeros(itrs, 1);
gradients   = cell(itrs, 1);

% Calculations for first permutation
le = getEmptyGraphLL(n, theta);
ge = getEmptyGraphGrad(n, theta);

% Calculate log likelihood
likelihood = le;
for e=1:size(E,1)
    eLL = edgeLL(n, theta, p(E(e,1)), p(E(e,2)));
    likelihood = likelihood - log(1 - exp(eLL)) + eLL;
end
% Calculate gradient
gradUpdate = ge;
for e=1:size(E,1)
    eGrad = edgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
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
            gradUpdate = ge;
            for e=1:size(E,1)
                eGrad = edgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
                gradUpdate = gradUpdate + eGrad;
                % gradUpdate = gradUpdate - log(1 - eGrad) + log(eGrad);
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

end
