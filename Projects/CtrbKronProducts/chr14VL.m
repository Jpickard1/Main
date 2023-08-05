%% Controlling Chromosome 14 with Kronecker Products
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu

clear; close all; clc

T = csvread("chr14_100kb_sparse_bin.csv");
T = T + 1;
S = spconvert(T);
A = full(S);

n = 9;
A = A(1:n^2,1:n^2);

[B,C] = nearestKroneckerProduct(A,[9 9], [9 9]);

clear S T

%% Plot factorization
figure;
tiledlayout(1,4);
nexttile;
imagesc(log10(A)); title('Chromosome 14 100kb');
nexttile;
imagesc(log10(B)); title('Kron Factor 1');
nexttile;
imagesc(log10(C)); title('Kron Factor 2');
nexttile;
imagesc(log10(kron(B,C))); title('Kron Rank 1 Approximation')

%% Svd factorization
[U,S,V] = svds(A,1)

figure;
tiledlayout(1,2);
nexttile;
imagesc(log10(A)); title('Chromosome 14 100kb');
nexttile;
imagesc(log10(S * U * V')); title('SVD Rank 1');

%% Controllability test

% Identify MCN on B and C

[CB, BiB] = greedyMCN(B);
[CC, BiC] = greedyMCN(C);

[B1] = inputMat(BiB, n);
[B2] = inputMat(BiC, n);

BB = kron(B1,B2)

CTRB = ctrb(A,BB);

CTRB(isnan(CTRB)) = 0;
CTRB(isinf(CTRB)) = 0;

rank(CTRB(:,1:100))

rank(ans)

