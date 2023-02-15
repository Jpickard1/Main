function [O, E] = observeSystemCE(HG,numGroups)
%OBSERVESYSTEM 
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023

[n, m] = size(HG.IM);

% Select hypergraph
popularity = hyperedgePopularity(HG.IM);
popS = sort(popularity);

idxs = round(1:m/(numGroups+1):m);

O = cell(numGroups, 1);
E = zeros(1, length(popS));
for g=1:numGroups
    EprimeBin = (popularity < popS(idxs(g))) + (popularity > popS(idxs(min(g + 1, end))));
    Eprime = find(EprimeBin ~= 0);
    IMp = zeros(n,m);
    IMp(:,Eprime) = HG.IM(:,Eprime);
    O{g} = Hypergraph('IM', IMp);
    E = E + EprimeBin;
end

end
