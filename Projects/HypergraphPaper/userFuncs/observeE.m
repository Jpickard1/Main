function [HGp] = observeE(HG, r, mp)
%OBSERVEE
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

IM = HG.IM;
[n, m] = size(IM);

% Bias and normalize popularity (line 2 in tex)
popularity = hyperedgePopularity(IM);
emptyEdges = find(popularity == 0);
popularity = popularity .^ r;
popularity(emptyEdges) = 0;

while isinf(sum(popularity))
    popularity = popularity / median(popularity);
end

% If there aren't enough hyperedges with nonzero weights to be selected,
% then the remainder are selected randomly.
%{
if length(find(popularity ~= 0)) < mp
    numRand = mp - length(find(popularity ~= 0));
    idxs = datasample(find(popularity == 0), numRand, 'Replace', false, 'Weights', ones(length(find(popularity == 0)), 1));
    popularity(idxs) = 1;
end
%}

% Normalize popularity
popularity = popularity / sum(popularity);

% Select hyperedges according to distribution (line 3 in tex)
Eprime = datasample(1:m,mp,'Replace',false,'Weights',popularity);

% Extract hypergraph (lines 4,5 in tex)
IMp = zeros(n,m);
IMp(:,Eprime) = IM(:,Eprime);
HGp = Hypergraph('IM', IMp);

end

%% Testing
%{
HG = HAT.uniformErdosRenyi(20,100,3)
IM = full(HG.IM);
r = 1
mp = 3
%}