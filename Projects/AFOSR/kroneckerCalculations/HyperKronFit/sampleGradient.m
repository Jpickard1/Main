function [likelihood, gradient]=sampleGradient(A, theta, debug)
% SAMPLEGRADIENT
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023
 
% TODO: rest
itrs = 50000; debug = false;

n0 = size(theta,1);
n = size(A,1);
kronExp = log(n) / log(n0);

% precompute edge set as a matrix
[e1, e2] = find(A ~= 0);
E = [e1 e2];

p = firstPermutation(A, theta); % p = 1:n; % For debugging purposes
if debug; disp(p); end
likelihoods = zeros(itrs, 1);
gradients   = cell(itrs, 1);

% Calculations for first permutation
le = getEmptyGraphLLapx(n, theta);
ge = getEmptyGraphGrad(n, theta);

% Calculate log likelihood
likelihood = le;
for e=1:size(E,1)
    eLL = edgeLLapx(n, theta, p(E(e,1)), p(E(e,2)));
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
    [pnew, accept, u, v] = nextPermutation(A, theta, p);

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

            % Get updated edges
            changedEdges = find(sum((E == u) + (E == v), 2) > 0);

            % Calculate log likelihood
            likelihood = likelihoods(t-1);
            for i=1:length(changedEdges)
                eLLold = edgeLLapx(n, theta, p(E(changedEdges(i),1)), p(E(changedEdges(i),2)));
                eLLnew = edgeLLapx(n, theta, pnew(E(changedEdges(i),1)), pnew(E(changedEdges(i),2)));
                likelihood = likelihood - ( -log(1 - exp(eLLold)) + ...
                    eLLold) + ( -log(1 - exp(eLLnew)) + eLLnew);
            end
            likelihoods(t) = likelihood;            
            gradUpdate = gradients{t-1};
            for i=1:length(changedEdges)
                eGradold = edgeGradient(n, theta, p(E(changedEdges(i),1)), p(E(changedEdges(i),2)));
                eGradnew = edgeGradient(n, theta, pnew(E(changedEdges(i),1)), pnew(E(changedEdges(i),2)));
                gradUpdate = gradUpdate - eGradold + eGradnew;
            end
            gradients{t} = real(gradUpdate);
            % likelihood = le;
            % for e=1:size(E,1)
            %     eLL = edgeLL(n, theta, p(E(e,1)), p(E(e,2)));
            %     likelihood = likelihood - log(1 - exp(eLL)) + eLL;
            % end
            % likelihoods(t) = likelihood;
        
            % Calculate gradient
            %   I could speed this loop up by summing only over the sets of
            %   edges that have changed (i.e. edges with endpoints u and v)
            % gradUpdate = ge;
            % for e=1:size(E,1)
            %     eGrad = edgeGradient(n, theta, p(E(e,1)), p(E(e,2)));
            %     gradUpdate = gradUpdate + eGrad;
            % end
            % gradients{t} = real(gradUpdate);

        % Save updated calculations
        likelihoods(t) = likelihood;
        gradients{t} = real(gradUpdate);
    else
        likelihoods(t) = likelihoods(t-1);
        gradients{t} = gradients{t-1};
    end
    p = pnew;
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
