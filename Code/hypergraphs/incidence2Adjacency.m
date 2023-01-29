function A = incidence2Adjacency(W)
%INCIDENCE 2 ADJACENCY 
%
%   This function transforms a k-way hypergraph incidence matrix to an
%   adjacanecy tensor.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 16, 2022

[n, e] = size(W);
k = max(sum(W));

% idxs = n * ones(k, 1);
A = symtensor(@zeros, k, n);

for i=1:e
    vxc = find(W(:,i) > 0);
    A(vxc) = W(vxc(1), i);
end

end

