%% Controllability of Kronecker Products
%
%   Aim: Mesbahi and Hao establish conditions where the controllability of
%   Kroncker product networks can be determined from the controllability of
%   the constituent systems according to a PBH style conditions; however,
%   niether work evaluated how often that condition is relevant in real
%   world systems. In this work, I perform the square Kronecker-rank one
%   approximation of several networks, determine the MCN on each
%   constituent system, and verify if the Kronecker product of the MCN
%   induces controllability on the composite system.


%% Chromosome 14 at 100kb Resolution

% Read file
T = csvread("chr14_100kb_sparse_bin.csv"); T = T + 1; S = spconvert(T);
% Form adjacency matrix
A = full(S); n = 9; A = A(1:n^2,1:n^2);
A = (A > median(A,'all')); % figure; imagesc(A)
A = double(A);
% A approx = B kron C
[B,C] = nearestKroneckerProduct(A,[9 9], [9 9]); clear S T

[CtrB, Bvec] = greedyMCN(B);
[CtrC, Cvec] = greedyMCN(C);

inB = inputMat(Bvec, n);
inC = inputMat(Cvec, n);

rank(ctrb(A, kron(inB,inC)))

[~,~,wB] = eig(B)
[~,~,wC] = eig(C)

%% Erdo Renyi Networks

n = 10;
A = randi([0 10], n^2); A = double(A == 1);

[B,C] = nearestKroneckerProduct(A,[n n], [n n]); clear S T

[CtrB, Bvec] = greedyMCN(B);
[CtrC, Cvec] = greedyMCN(C);

inB = inputMat(Bvec, n);
inC = inputMat(Cvec, n);

rank(ctrb(A, kron(inB,inC)))