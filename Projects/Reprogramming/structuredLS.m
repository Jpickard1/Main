function [C] = structuredLS(A,B,S)
%STRUCTUREDLS Least Squares with Structured Sparsity
%
%   This function solves the minimization problem for C:
%
%       min ||A-B*C|| such that (C ~= 0) == S
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 14, 2023

C = zeros(size(S));
for i=1:size(C,2)
    N = find(S(:,i) ~= 0);
    Ci = B(:,N) \ A(:,i);
    C(N,i) = Ci;
end


end

