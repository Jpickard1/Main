function [likelihood, gradient]=sampleGradient(A, theta, itrs, firstPermItrs, debug, directed, E)
% SAMPLEGRADIENT
% 
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

k = ndims(A);
n0 = size(theta,1);
n = size(A,1);
kronExp = log(n) / log(n0);

% precompute edge set as a matrix
E = getEdgesFromAdj(A);

p = firstPermutation(A, theta, firstPermItrs, debug, E);

% p = 1:n; % For debugging purposes
if debug; disp(p); end % TODO: this should probably be remved at least once I verify the code can work

likelihoods = zeros(itrs, 1);
gradients   = cell(itrs, 1);

% Calculations for first permutation
le = getEmptyHypergraphLLapx(n, theta);
ge = getEmptyHypergraphGrad(n, theta, directed); 

likelihood = le;
gradUpdate = ge;
for e=1:size(E,1)
    % Calculate log likelihood
    eLL = hedgeLLapx(n, theta, p(E(e,:)));
    likelihood = likelihood - log(1 - exp(eLL)) + eLL;
    % Calculate gradient
    eGrad = hedgeGradient(n, theta, p(E(e,:)), eLL, directed);
    gradUpdate = gradUpdate + eGrad;
end
likelihoods(1) = likelihood;
gradients{1} = gradUpdate;

accepted = 0;
for t=2:itrs
    % if debug; disp(p); end
    % Generate sample permutation
    [pnew, accept, u, v] = nextPermutation(A, theta, p, E);

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
            % This is not so fast, its actually probably a pretty slow line
            changedEdges = find(sum((E == u) + (E == v), 2) > 0);

            % Calculate log likelihood
            likelihood = likelihoods(t-1);
            gradUpdate = gradients{t-1};
            for i=1:length(changedEdges)
                % log likelihoods
                eLLold = hedgeLLapx(n, theta, p(E(changedEdges(i),:)));
                eLLnew = hedgeLLapx(n, theta, pnew(E(changedEdges(i),:)));
                likelihood = likelihood - ( -log(1 - exp(eLLold)) + ...
                    eLLold) + ( -log(1 - exp(eLLnew)) + eLLnew);
                % gradients
                eGradold = hedgeGradient(n, theta, p(E(changedEdges(i),:)), eLLold, directed);
                eGradnew = hedgeGradient(n, theta, pnew(E(changedEdges(i),:)), eLLnew, directed);
                gradUpdate = gradUpdate - eGradold + eGradnew;
            end
            likelihoods(t) = likelihood;
            gradients{t} = real(gradUpdate);
            % for i=1:length(changedEdges)
            %     eGradold = hedgeGradient(n, theta, p(E(changedEdges(i),:)), directed);
            %     eGradnew = hedgeGradient(n, theta, pnew(E(changedEdges(i),:)), directed);
            %     gradUpdate = gradUpdate - eGradold + eGradnew;
            % end
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
