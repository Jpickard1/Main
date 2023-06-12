function [gradient] = hedgeGradient(n, theta, idxs)
%HEDGEGRADIENT Evaluates the gradient of theta with respect to a hyperedge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

theta = theta / sum(theta, "all");

n0 = size(theta,1);
kronExp = log(n) / log(n0);

eLL = hedgeLL(n, theta, idxs);
noEdgeLL = log(1 - exp(eLL));

% Count the number of times an entry of theta is used
count = zeros(size(theta));
for i=1:kronExp
    idx = mod(floor((idxs - 1) / n0^(i - 1)), n0) + 1;
    idx = num2cell(idx);
    count(idx{:}) = count(idx{:}) + 1;
end

% Calculate gradient of log likelihood function
gradient = zeros(size(theta));
E = getAllHyperedges(n0, k);
for i=1:size(E,1)
    idx = E(i,:); idx = num2cell(idx);

    % Count the number of times (i,j) was used
    c = count(idx{:});

    % Set gradient if edge is present
    posGrad = c / exp(theta(idx{:}));

    % Set gradient if edge is not present
    negGrad = getNoHedgeDLL2(theta, count, idx, eLL);

    % Calculate gradient update
    gradUpdate = posGrad - negGrad;
    
    % Update all symmetric entries of the gradient tensor
    idxp = perms(cell2mat(idx));
    idxp = unique(idxp, 'rows');
    for j=1:size(idxp,1)
        pidx = num2cell(idxp(j,:));
        gradient(pidx{:}) = gradUpdate;
    end
end

end