function [likelihood, gradient, p]=sampleGradient(theta, itrs, firstPermItrs, E, p, pFixed)
% SAMPLEGRADIENT
% 
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

k = size(E,2);
n0 = size(theta,1);
n = max(max(E));
kronExp = log(n) / log(n0);

if nargin == 4
    p = randperm(n);
end
if nargin < 6
    pFixed = false;
end

if ~pFixed
    p = genPermutation(p, theta, firstPermItrs, E);
end

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
    if pFixed
        break
    end
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

if pFixed
    likelihood = likelihoods(1);
    gradient = gradients{1};
else
    likelihood = mean(likelihoods);
    gradient = zeros(size(theta));
    for i=1:itrs
        gradient = gradient + gradients{i};
    end
    gradient = gradient ./ itrs;
    fprintf("Accepted: %d / %d \n", accepted, itrs);
end

end
