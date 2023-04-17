function [S] = shuffleMatrix(p, q)
%SHUFFLEMATRIX This makes the perfect shuffle matrix that can be used to
% make the Kronecker product associative.
%
%EXAMPLE
%
%   A = rand(3,3);
%   B = rand(4,4);
%   S = shuffleMatrix(3, 4)
%
%   kron(A, B) == S' * kron(B, A) * S
%
%REFERENCE
%   1. Kronecker Products and Shuffle Algera, Marc Davio (1981)
%   2. https://en.wikipedia.org/wiki/Kronecker_product
%
% Auth: Joshua B. Pickard
%       jpic@umich.edu
% Date: April 8, 2023

I = eye(p*q,p*q);

S = I(1:q:p*q,:);
for i = 2:q
    S = [S; I(i:q:p*q,:)];
end

end

