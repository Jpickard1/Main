function [B] = TenKronFit(A, tol)
%TENKRONFIT This function performs a tensor Kronecker decomposition of an
%   arbitrary size tensor A such that A = superkron(B,B)
%
% Algorithms:
%   1. Initialize B
%   2. while not converged
%       3. evaluate gradient
%       4. update B
%       5. check convergence
%   6. return B
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 18, 2023

clear; close all; clc;
A = rand(4,4,4);

s = size(A);
B = rand(s / 2);

converged = false;
while ~converged
    G = getGradient()
    converged = checkConverged(A, B);
end

end

function C = checkConverged(A, B, tol)
    BB = superkron(B, B);
    n = norm(A-BB, 'fro');
    C = false;
    if n < tol; C = true; end
end

function [G,  L] = getGradient(B, A)
    T = 10; % TODO: Change this
    grad = zeros(T,1);
    likelihood = zeros(T,1);
    for i=1:T
        sigma = samplePermutation(A, B);
        likelihood(i) = 
        grad(i) = 
    end
    G = mean(grad);
    L = mean(lieklihood);
end
