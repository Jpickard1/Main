%% Is Observable Kronecker Hypergraphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 4, 2023

clear; clc; close all;

V = 5;
k = 3;
O = [1 2 3];
HG = hyperstar(V, k);

C = zeros(length(O), V);
for i=1:length(O); C(i,i) = 1; end
% C

A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
n = size(A, 1);             % Get n
k = length(size(A));

% Unfold A
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

% Compute J vectors (equations 5 and 6 in overleaf document)
% This is the time computationally expensive part of the program
J = cell(n,1);
J{1} = sym('x_%d',[n 1]);       % J0 in latex
sVec = loadSymVecs(n, size(Amat,2));
J{2} = Amat * sVec;             % J1 in latex
P = sparse(Amat);               % Tracks terms AB1B2B... in equation 6
for i=2:n
    % getBp is equation 7 in overleaf document 
    P = P * getBp(sparse(Amat), i, k);
    sVec = loadSymVecs(n, size(P,2));
    J{i} = P * sVec;            % Equation 6
end

% Compute observability matrices for all vertices
symVars = symvar(sym('x_%d',[n 1]));  % Get symbolic variables 
Oi = cell(n,1);             % Compute first equality in equation 10
for i=1:n
    Oi{i} = C * J{i};      % Compute first equality in equation 10
end
Oimat = sym([n,n]);         
for i=1:n                   % Set symbolic matrix to save equation 10 for vertex vx
    for j=1:n               % Compute second equality in equation 10
        Oimat(i,j) = gradient(Oi{i}, symVars(j));
    end
end
Oimat;              % Save observability matrix for specific vertex

