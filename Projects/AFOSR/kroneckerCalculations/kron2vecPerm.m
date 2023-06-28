function [Ap] = kron2vecPerm(A)
%KRON2VECPERM Transforms a matrix from a kronecker equation to a matrix
%   vector equation i.e. if we have:
%
%       A = B kron C
%
%   then this function is used so that
%
%       kron2vecPerm(A) = B(:) * C(:)
%
%   currently this assumes A is square and B and C are equal sizes
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 28, 2023

if isnumeric(A)
    Ap = zeros(size(A));
else
    Ap = sym('x', size(A));
end

n2 = size(A,1);
n = round(sqrt(n2));
for i=1:n
    for j=1:n
        AA = A( ((i-1)*n+1):i*n, ((j-1)*n+1:j*n))';
        Ap((i-1) * n + j,:) = AA(:); % reshape(A( ((i-1)*n+1):i*n, ((j-1)*n+1:j*n)), [1 n2]);
    end
end

end
