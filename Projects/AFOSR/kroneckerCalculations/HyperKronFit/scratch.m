%% Scratch
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 11, 2023

%% Debugging PPRtest2

n = 4;
A = rand(n,n);
theta = [0.9 0.8;
         0.6 0.1];

p=1:n;
u = 1; v = 2;
PPRtest2(p, theta, A, u, v)



p1lluC = 0;
for i=1:n
    ell = edgeLLapx(n, theta, i, 1);
    p1lluC = p1lluC + log(1 - exp(ell));
end
disp(p1lluC)

P = kron(theta, kron(theta, theta));
s = 0;
for i=1:n
    s = s - P(1,i) - 0.5 * (P(1,i) ^ 2);
end
disp(s)

s = 0;
for i=1:n
    s = s - P(i,1) - 0.5 * (P(i,1) ^ 2);
end
disp(s)

%%
clear; clc; close all;
n0 = 2; n = 4;
kronExp = log(n) / log(n0);
theta = rand(n0, n0); P = theta;
for i=2:kronExp
    P = kron(theta,P);
end

sum(sum(P))

sum(sum(theta)) ^ kronExp

sum(P(4,:))
sum(theta(2,:)) ^ kronExp

sum(P(2,:))
sum(theta(1,:)) * sum(theta(2,:))
