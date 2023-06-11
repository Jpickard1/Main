%% Kronecker Fit and Kronecker Generator Demo File
%
%   This file demonstrats the funcitonality of differnet components
%   required to work with Kronecker graphs and hypergraphs.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 8, 2023

clear; clc; close all;

%% Graph Alignment
kronExp = 3;
eps = 0.1;
theta = [1 1 1;
         0 1 0
         1 0 1];
theta = theta - eps; theta(theta < 0) = eps;

P = theta;
for i=1:kronExp
    P = kron(theta,P);
end
A = (P > 0.25);
tic;
[p, pp] = firstPermutation(A, theta, 500000, false);
t1 = toc;

% Check graph alignment of A into P using the perms
n = size(A,1);
Ap = zeros(n,n);
for i=1:n
    for j=1:n
        Ap(i,j) = A(p(i), p(j));
    end
end

figure; 
subplot(1,3,1); imagesc(A); title('Kronecker Graph');
subplot(1,3,2); imagesc(P); title('Kronecker Expansion');
subplot(1,3,3); imagesc(Ap); title('Aligned Graph');

h = figure;
% h.Visible = 'off';
M(size(pp,1)) = struct('cdata',[],'colormap',[]);
for t = 1:size(pp,1)
    Ap = zeros(n,n);
    for i=1:n
        for j=1:n
            Ap(i,j) = A(pp(t,i), pp(t,j));
        end
    end
    imagesc(Ap);
    drawnow
    M(t) = getframe;
end
% h.Visible = 'on';
% figure; movie(M); title('Graph Alignment Process');

%% Kronecker Graph Generation

theta = [1 1 1
         1 1 0
         1 0 1];
eps = 0.01;
theta = theta - eps; theta(theta < 0) = eps;
n0 = size(theta, 1);
kronExp = 4;
numE = 50 * n0^kronExp;

theta = theta / sum(sum(theta));
E = kronGen(theta, kronExp, numE);
n = n0^kronExp;
A = sparse(n, n);
for i=1:size(E,1)
    A(E(i,1), E(i,2)) = 1;
end
A = full(A);
figure; imagesc(A)


%% KronFit from KronGen

theta = rand(3,3);
eps = 0.05;
theta = theta - eps; theta(theta < 0) = eps;
n0 = size(theta, 1);
kronExp = 5;
numE = 200 * n0^kronExp;

theta = theta / sum(sum(theta));
E = kronGen(theta, kronExp, numE);
n = n0^kronExp;
A = sparse(n, n);
for i=1:size(E,1)
    A(E(i,1), E(i,2)) = 1;
end
A = full(A);
figure; imagesc(A)

%% Using A from above
[thetaLearned, likelihoods] = NaiveKronFit(A, true, true, 3);
[thetaLearned2, likelihoods2] = NaiveKronFit(A, true, true, 3, thetaLearned);
[thetaLearned3, likelihoods3] = NaiveKronFit(A, true, true, 3, thetaLearned2);

ll = [likelihoods; likelihoods2; likelihoods3];
figure; plot(ll);

thetaLearned3 = thetaLearned3 / sum(sum(thetaLearned3));
E = kronGen(thetaLearned3, kronExp, numE);
n = n0^kronExp;
A = sparse(n, n);
for i=1:size(E,1)
    A(E(i,1), E(i,2)) = 1;
end
A = full(A);
figure; imagesc(A)

%%
theta = [0.1385    0.1385    0.1385
         0.1385    0.1385    0.0154
         0.1385    0.0154    0.1385]


