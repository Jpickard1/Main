%% Kronecker Fit and Kronecker Generator Demo File
%
%   This file demonstrats the funcitonality of differnet components
%   required to work with Kronecker graphs and hypergraphs.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 8, 2023

clear; clc; close all;

%% Kronecker Graph Generation, Fitting, and Regeneration

% Set parameters
theta0 = [0.9 0.6;
          0.6 0.1];
n0 = 2;
kronExp = 14;
numE = 5 * n0^kronExp;
theta00 = theta0 / sum(theta0, 'all');

% Generate kronecker graph
E = kronGen(theta00, kronExp, numE);
n = n0^kronExp;
A0 = sparse(n, n);
for i=1:size(E,1)
    A0(E(i,1), E(i,2)) = 1;
end
A0 = full(A0);

% Permute Kronecker Graph
E = zeros(sum(A0, 'all'), 2);
[E(:,1), E(:,2)] = find(A0 > 0);
p = randperm(n);
Ap = sparse(n,n);
for e=1:size(E,1)
    Ap(p(E(e,1)), p(E(e,2))) = 1;
end

% Learn parameters from permuted graph structure
thetaG = rand(2,2);
[thetaLearned, likelihoods] = NaiveKronFit(Ap, true, true, 2, thetaG, 20);

% Relearned kronecker graph
thetaLearned0 = thetaLearned / sum(thetaLearned, 'all');
E = kronGen(thetaLearned0, kronExp, numE);
n = n0^kronExp;
Ar = sparse(n, n);
for i=1:size(E,1)
    Ar(E(i,1), E(i,2)) = 1;
end
Ar = full(Ar);

figure; nplots = 6;
subplot(1,nplots,1); imagesc(theta0); title('Initiator Parameters');
subplot(1,nplots,2); imagesc(A0); title('Kronecker Graph');
subplot(1,nplots,3); imagesc(Ap); title('Permuted Kronecker Graph');
subplot(1,nplots,4); imagesc(thetaG); title('Initial Guess of Parameters');
subplot(1,nplots,5); imagesc(thetaLearned); title('Estimated Parameters');
subplot(1,nplots,6); imagesc(Ar); title('Regenerated Kronecker Graph');




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
[p, pp] = firstPermutation(A, theta, 175000, false);
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

theta = [1 0 1
         0 1 0
         1 0 1];
eps = 0.01;
theta = theta - eps; theta(theta < 0) = eps;
n0 = size(theta, 1);
kronExp = 5;
numE = 50 * n0^kronExp;

theta = theta / sum(sum(theta));
E = kronGen(theta00, kronExp, numE);
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
theta0 = rand(3,3);
[thetaLearned, likelihoods] = NaiveKronFit(A, true, true, 3, theta0, 100);
[thetaLearned2, likelihoods2] = NaiveKronFit(A, true, true, 3, thetaLearned, 25);
[thetaLearned3, likelihoods3] = NaiveKronFit(A, true, true, 3, thetaLearned2, 25);
[thetaLearned4, likelihoods4] = NaiveKronFit(A, true, true, 3, thetaLearned3, 25);
[thetaLearned5, likelihoods5] = NaiveKronFit(A, true, true, 3, thetaLearned4, 100);
[thetaLearned6, likelihoods6] = NaiveKronFit(A, true, true, 3, thetaLearned5, 200);
[thetaLearned7, likelihoods7] = NaiveKronFit(A, true, true, 3, thetaLearned6, 50);

ll = [likelihoods; likelihoods2; likelihoods3; likelihoods4; likelihoods5; likelihoods6; likelihoods7];
figure; plot(ll);

thetaLearnedF = thetaLearned7 / sum(sum(thetaLearned7));
E = kronGen(theta, kronExp, numE);
n = n0^kronExp;
A = sparse(n, n);
for i=1:size(E,1)
    A(E(i,1), E(i,2)) = 1;
end
A = full(A);
figure; imagesc(A)


thetaLearnedF = thetaLearned7 / sum(sum(thetaLearned7));
E = kronGen(thetaLearnedF, kronExp, numE);
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


