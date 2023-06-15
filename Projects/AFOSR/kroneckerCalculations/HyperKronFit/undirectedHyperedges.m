function [EE] = undirectedHyperedges(E)
%UNDIRECTEDHYPEREDGES
%
%   E is a m x k matrix for a k-uniform hypergraph on m vertices. Each
%   hyperedge is given in a particular order. The purpose of this function
%   is to return all possible orders of vertices in each hyperedge.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 15, 2023

E = unique(sort(E, 2), 'rows');
k = size(E, 2);
rpts = factorial(k);
EE = zeros(size(E,1) * rpts, k);

for e=1:size(E,1)
    hedge = E(e,:);
    EE(((e-1) * rpts + 1):(e * rpts),:) = perms(hedge);
end

EE = unique(EE, 'rows');

end