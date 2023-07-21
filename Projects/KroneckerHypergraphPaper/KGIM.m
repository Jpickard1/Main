%% Kronecker Graphs Incidence Matrices
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 14, 2023


clear; close all; clc
n = 10;
e = 15;
A1 = erdos_renyi_network(n, e); A1 = double(A1);
A2 = erdos_renyi_network(n, e); A2 = double(A2);

[x, y] = find(A1 == 1); E1 = [x y];
[x, y] = find(A2 == 1); E2 = [x y];

I1 = E2IM(E1, n);
I2 = E2IM(E2, n);

A = kron(A1,A2);
I = kron(I1,I2);

sum(A == I * I' - diag(sum(I,2)), 'all')

sum(A1 == I1 * I1' - diag(sum(A1)),'all')

%%
clear all; clc; close all;
n = 2;
B = rand(n,n,n);
C = rand(n,n,n);
A = superkron(B,C);

Ar = reshape(A,[n^2 n^4]);
BrCr = kron(reshape(B,[n n^2]), reshape(C,[n n^2]));

Ar == BrCr

%% 
clear all; close all; clc

for i=1:1000
    disp(i)
    A = rand(2,2,2);
    ea = heig(A);
    B = rand(2,2,2);
    eb = heig(B);
    C = superkron(A,B);
    ec = heig(C);
    
    if min(ea) > 0 && min(eb) > 0 && min(ec) < 0
        disp(A);
        disp(B);
        disp(C);
        break;
    end
end

%%

clear
n = 
A = rand(n,n,n);
tic;
ea = heig(A);
ta = toc
disp(ta);
B = rand(n,n,n);
tic;
eb = heig(B);
tb = toc
disp(tb);
C = superkron(A,B);
tic;
ec = heig(C);
tc = toc
disp(tc)

