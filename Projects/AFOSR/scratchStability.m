%% Stability of Kronecker Continuous and Discrete Linear Systems
%
%   I am interested in the following 2 equations:
%
%       (1) dx/dt  = Ax
%       (2) x(t+1) = Ax
%
%   where A is an n dimensional matrix and x is an n dimensional vector.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 4, 2023

clear; clc; close all;

T = 100; n = 2;
A = rand(n,n);
ceq = @(t,x) -A * x;
X0 = rand(n,1);

[t, X1] = ode45(ceq,[1:0.01:T], X0);
figure; plot(t, X1); hold on;

% As a difference equation
At = exp(A * 0.01);
X2 = zeros(n,T); X2(:,1) = X0;
for t=2:size(X2,2)
    X2(:,t) = At * X2(:,t-1) + X2(:,t-1);
end
scatter(t, X2(1,:));
scatter(t, X2(2,:));

X1 = zeros(n,T);
X2 = zeros(n,T);
X1(:,1) = rand(1,n);

for t=2:T
    X1(:,t) = A * X1(:,t-1) + X1(:,t-1);
end

figure; plot(X')


A

