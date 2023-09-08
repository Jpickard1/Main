function [G] = observabilityGrammarian(A, C, maxN)
%OBSERVABILITYGRAMMARIAN
%
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

if nargin == 2
    maxN = 10;
end

CC = C*C';

G = zeros(size(A));
for i=0:maxN
    Ai = A^i;
    G = G + Ai' * CC * Ai;
end

end

