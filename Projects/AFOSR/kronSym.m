function [K] = kronSym(A, B)
%KRONSYM Kronecker Producet for symbolic matrices
%
%   I think this code is unnecessary. I wrote this because I was having
%   trouble doing the kronecker product of 2 symbolic matrices where the
%   matrices were tensor() toolbox objects, but it seems that the trouble
%   was really how the tensor objects are handled. 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 5, 2023

s1 = size(A);
s2 = size(B);

n1 = s1(1);
m1 = s1(2);
n2 = s2(1);
m2 = s2(2);

%%
K = sym('x', [n1 * n2, m1 * m2]);
for i=1:n1
    for j=1:m1
        K((i-1)*n2+1:i*n2, (j-1)*m2+1:j*m2) = A(i,j) * B;
    end
end

end

%{
    if isSymType(A(1,1), "expression") || isSymType(A(1,1), "variable")
    else
        K = zeros
    end
%}