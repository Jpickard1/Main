function [E] = getHyperedges(A)
%GETHYPEREDGES 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

k = ndims(A);
idxs = cell(1, k);
linIdx = find(A > 0);
[idxs{:}] = ind2sub(size(A), linIdx);
E = cell2mat(idxs);

end

