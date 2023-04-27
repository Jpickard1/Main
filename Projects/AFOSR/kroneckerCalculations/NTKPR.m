function [K,N] = NTKPR(A, R)
%NTKPR
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 18, 2023

N = zeros(R+1,1);
N(1) = norm(A,'fro');
disp(norm(A, 'fro'));
K = cell(R,2);
for r=1:R
    [K{r,1}, K{r,2}] = NTKP(A);
    A = A - superkron(K{r,1}, K{r,2});
    disp(norm(A, 'fro'));
    N(r+1) = norm(A,'fro');
end

end

