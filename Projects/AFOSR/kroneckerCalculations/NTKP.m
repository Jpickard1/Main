function [B, C] = NTKP(A)
%NTKP Nearest Tensor Kronecker Product
%
% This function returns tensors B and C that solve the Nearest Tensor
% Kronecker Product problem descibed by the minimization:
%
%                          min |A - kron(B, C)|
%
% The algorithm has 4 steps:
%   1. Ak       = Unfold A
%   2. Ap       = Permute columns of Ak
%   3. [Bm, Cm] = Nearest Kronecker Product of Ap
%   4. [B, C]   = reshape(Bm), reshape(Cm)
%
% Auth: Joshua B. Pickard
%       jpic@umich.edu
% Date: April 9, 2023

%% Start
s = size(A);
if length(unique(s)) ~= 1
    error('A is not k-uniform');
end
n1 = sqrt(s(1));
k = size(s, 2);

%%   1. Ak       = Unfold A
Ak = reshape(A, [size(A,1), numel(A) / size(A,1)]);

%%   2. Ap       = Permute columns of Ak
n = 1:n1;
p = [];
for i=1:n1
    p = [p ((n1^2)*(i-1) + n)];
end
rp = [];
for i=1:n1
    rp = [rp (n1*(i-1) + p)];
end
permutation = zeros(size(Ak, 2), 1);
for i=1:n1
    permutation(1+(i-1)*length(rp):(i)*length(rp)) = (n1^k)*(i-1) + rp;
%    permutations = [permutations ((n1^k)*(i-1) + rp)];
end
Ap = Ak(:, permutation);

%   3. [Bm, Cm] = Nearest Kronecker Product of Ap
sq = sqrt(s);
sB = [sq(1) (prod(sq)/sq(1))];    % size of B matrix
sC = [sq(1) (prod(sq)/sq(1))];    % size of B matrix
[Bv, Cv] = nearestKroneckerProduct(Ap, sB, sC);

%   4. [B, C]   = reshape(Bm), reshape(Cm)
B = reshape(Bv, sqrt(s));
C = reshape(Cv, sqrt(s));

end

