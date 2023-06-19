function [likelihood, gradient]=sampleGradient(theta, itrs, firstPermItrs, E)
% SAMPLEGRADIENT
% 
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

k = size(E,2);
n0 = size(theta,1);
n = max(max(E));
kronExp = log(n) / log(n0);

% precompute edge set as a matrix
p = randperm(n);
p = genPermutation(p, theta, firstPermItrs, E);

likelihoods = zeros(itrs, 1);
gradients   = cell(itrs, 1);

% Calculations for first permutation
le = emptyLL(n, theta);
ge = emptyDLL(n, theta); 

likelihood = le;
gradUpdate = ge;
for e=1:size(E,1)
    % Calculate log likelihood
    eLL = hedgeLL(n, theta, p(E(e,:)));
    likelihood = likelihood - log(1 - exp(eLL)) + eLL;

    % Calculate gradient
    eGradPos = hedgeDLL(n, theta, p(E(e,:)));
    eGradNeg = noHEdgeDLL(n, theta, p(E(e,:)));
    gradUpdate = gradUpdate + eGradPos - eGradNeg;
end
likelihoods(1) = likelihood;
gradients{1} = gradUpdate;

accepted = 0;
for t=2:itrs
    % Generate sample permutation
    [pnew, accept, idxs] = genPermutation(p, theta, 1, E);

    % Update 
    if accept
        % Update number of accepted permutations
        accepted = accepted + 1;

        % Get updated edges
        % This is not so fast, its actually probably a pretty slow line
        eIdxs = find(sum(ismember(E,idxs), 2));
        
        % Calculate log likelihood
        likelihood = likelihoods(t-1);
        gradUpdate = gradients{t-1};
        for i=1:length(eIdxs)
            % log likelihoods
            eLLold = hedgeLL(n, theta, p(E(eIdxs(i),:)));
            eLLnew = hedgeLL(n, theta, pnew(E(eIdxs(i),:)));
            likelihood = likelihood - ( -log(1 - exp(eLLold)) + eLLold) + ( -log(1 - exp(eLLnew)) + eLLnew);

            % gradients
            eDLLold = hedgeDLL(n, theta, p(E(eIdxs(i),:))) - noHEdgeDLL(n, theta, p(E(eIdxs(i),:)));
            eDLLnew = hedgeDLL(n, theta, pnew(E(eIdxs(i),:))) - noHEdgeDLL(n, theta, pnew(E(eIdxs(i),:)));
            gradUpdate = gradUpdate - eDLLold + eDLLnew;
        end
        % Save updated calculations
        likelihoods(t) = likelihood;
        gradients{t} = real(gradUpdate);
    else
        likelihoods(t) = likelihoods(t-1);
        gradients{t} = gradients{t-1};
    end
    p = pnew;
end

likelihood = mean(likelihoods);
gradient = zeros(size(theta));
for i=1:itrs
    gradient = gradient + gradients{i};
end
gradient = gradient ./ itrs;

fprintf("Accepted: %d / %d \n", accepted, itrs);
% fprintf("CurrentLL: %d\n", l);
% fprintf("Gradient Updates:\n");
% for i=1:numel(theta)
%     fprintf("    %d]  %f \t Grad:  %f\n", i, theta(i), gradient(i));
% end
% fprintf(' \n');

end
