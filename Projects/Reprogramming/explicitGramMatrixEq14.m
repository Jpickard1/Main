function [G] = explicitGramMatrixEq14(A,x0,maxI)
%GRAMMATRIXEQ14 Computes the Gram Matrix according to Equation 14
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 12, 2023

if nargin == 2
    maxI = 10;
end

G = zeros(numel(x0), numel(x0));

Ai = eye(size(A));
for i=0:maxI
    G = G + (Ai * x0) * (x0' * Ai');
    Ai = Ai * A;
end

end