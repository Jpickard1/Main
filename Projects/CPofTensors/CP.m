function [] = CP(T)
%CP constructs the characteristic polynomial of tensor T
%
%   Tx^{k-1} <==> Tpx^{[k-1]}
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2023
T = A
Tp = reshape(T, [size(T,1) (numel(T)/size(T,1))]);
X = sym('xv_%d', [size(T,1), 1])
Tp * kron(X,X)

end

