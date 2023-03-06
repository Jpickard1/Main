function [observable] = kroneckerObservable(A1,A2,B1,B2)
%KRONECKEROBSERVABLE This function tests if a graph A is observable
%   (controllable) with vertices B given kronecker gractors where
%
% >>    A = kron(A1, A2)
% >>    B = kron(B1, B2)
%
%   according to Theorem 6 in Kronecker Procut of Networked Systems and
%   their Approximates by Airlie Chapman and Mehran Mesbahi
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 6, 2023

n1 = size(A1);
n2 = size(A2);

% This function uses controllability tests, but it is the same as
% observability tests; however, this does require using the transposed
% matrix.
B1 = B1';
B2 = B2';   

r1 = rank(ctrb(A1,B1));
r2 = rank(ctrb(A2,B2));
if any([r1 ~= n1, r2 ~= n2])
    observable = false;
    return;
end

% Construct kronecker graph
A = kron(A1, A2);

% Diagonalize
[v, d] = eig(A);
isD = isdiag(inv(v) * A * v);

% Check if left eigenvectors of Ai are orthognal to Bi
[~,~,W1] = eig(A1);
[~,~,W2] = eig(A2);
s1 = W1' * B1;
s2 = W2' * B2;

if all([isD, sum(find(s1 == 0)) == 0, sum(find(s2 == 0)) == 0])
    observable = true;
else
    observable = false;
end
