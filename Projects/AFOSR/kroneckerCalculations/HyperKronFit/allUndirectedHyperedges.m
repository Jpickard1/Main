function [E] = allUndirectedHyperedges(n, k)
%GETALLHYPEREDGES Returns the indices of all undirected hyperedges
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

linIdx = 1:n^k;                         % Get all linear indices
idxs = cell(1, k);                      % Linear to subindicies
[idxs{:}] = ind2sub(n * ones(1, k), linIdx');
E = cell2mat(idxs);
E = sort(E')';                          % Sort the elements of each hyperedge
E = unique(E, 'rows');                  % Remove duplicate rows

end
