function [C] = kronSum(A, B)
%KRONSUM This function computes the Kronecker Sum of matrices A and B
% according to the equation:
%
%   C = kron(A, Ib) + kron(Ia,B)
%
% where Ib and Ia are identity matrices of appropriate sizes to ensure that
% the matrix addition is appropriate.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 4, 2023

na = size(A,1);
nb = size(B,1);

C = kron(A, eye(nb)) + kron(eye(na), B);

end

