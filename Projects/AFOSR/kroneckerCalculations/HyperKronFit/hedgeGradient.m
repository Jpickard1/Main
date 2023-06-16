function [gradient] = hedgeGradient(n, theta, idxs, directed)
%HEDGEGRADIENT Evaluates the gradient of theta with respect to a hyperedge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

% theta = theta / sum(theta, "all");

k = ndims(theta);
n0 = size(theta,1);
kronExp = log(n) / log(n0);

eLL = hedgeLLapx(n, theta, idxs);
% noEdgeLL = log(1 - exp(eLL));

% Count the number of times an entry of theta is used
count = zeros(size(theta));
for i=1:kronExp
    idx = mod(floor((idxs - 1) / n0^(i - 1)), n0) + 1;
    idx = num2cell(idx);
    count(idx{:}) = count(idx{:}) + 1;
end

% Calculate gradient of log likelihood function with respect to each
% parameter in theta
gradient = zeros(size(theta));
if ~directed     % Get list of entries in theta that need to be updated
    E = allUndirectedHyperedges(n0, k);
else
    E = allDirectedHyperedges(n0, k);
end
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
    if ~directed
        idxp = perms(cell2mat(idx));
        idxp = unique(idxp, 'rows');
        for j=1:size(idxp,1)
            pidx = num2cell(idxp(j,:));
            gradient(pidx{:}) = gradUpdate;
        end
    else
        gradient(idx{:}) = gradUpdate;
    end
end

end