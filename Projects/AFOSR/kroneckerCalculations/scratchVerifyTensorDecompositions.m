%% Verify the tensor decomposition results from the paper
%
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 26, 2023

clear; close all; clc;

n = 2;
m = 3;
B = rand(n,n,n);
C = rand(m,m,m);
A = superkron(B, C);

At = tensor(A);
Bt = tensor(B);
Ct = tensor(C);

acp = cp_als(tensor(A), 1);
bcp = cp_als(tensor(B), 1);
ccp = cp_als(tensor(C), 1);

