%% Simultaneous Kronecker Single Value Decomposition
%
%   This file tests the ability to perform simultaneous Kronecker SVD on 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 29, 2023

n = 4;                      % Set matrix size
n2 = n^2;
% Q = orth(rand(n));        % make orthonormal matrix
v = rand(n2,1);
v = v / norm(v);
Q = eye(n2) - 2 * v * v';    % Householder reflection

A = Q*diag(rand(n2,1))*Q';    % Make 2 similar matrices A and B
B = Q*diag(rand(n2,1))*Q';

[Ua, Sa, Va] = svd(A)
[Ub, Sb, Vb] = svd(B)

[Ba, Ca] = nearestKroneckerProduct(A, [(n) (n)], [(n) (n)])
[Bb, Cb] = nearestKroneckerProduct(B, [(n) (n)], [(n) (n)])

