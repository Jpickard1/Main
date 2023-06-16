function emptyGradient = getEmptyHypergraphGrad(n, theta, directed)
%GETEMPTYGRAPHGRAD
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

% theta = theta / sum(theta, 'all');

emptyGradient = zeros(size(theta));

n0 = size(theta,1);
k = length(size(theta));

kronExp = log(n) / log(n0);

s  = sum(theta, 'all');
sq = sum(theta .^ 2, 'all');

if ~directed
    E = allUndirectedHyperedges(n0, k);
else
    E = allDirectedHyperedges(n0, k);
end
for i=1:size(E,1)
    idx = E(i,:); idx = num2cell(idx);
    gradUpdate = -kronExp * (s ^ (kronExp - 1)) - kronExp * (sq ^ (kronExp - 1)) * theta(idx{:});
    if ~directed
        idxp = perms(cell2mat(idx));
        idxp = unique(idxp, 'rows');
        for j=1:size(idxp,1)
            pidx = num2cell(idxp(j,:));
            emptyGradient(pidx{:}) = gradUpdate;
        end
    else
        emptyGradient(idx{:}) = gradUpdate;        
    end
end

end