%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SYSTEM INTEGRATION                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Kronecker parameters (theta) and the graph alignment permutation (p)
%   are the only parameters in the system. To test the system integration,
%   we:
%       1) generate a kronecker hypergraph with initial theta and p
%       2) freeze theta (p)
%       3) estimate p (theta)
%       4) compare the estimated value of p (theta) to the true value used
%          to generate the system
%
% GRAPH ALIGNMENT 1
%   Purpose: Test code's ability to identify a correct mapping between 2
%            graphs using metropolis sampling algorithm
%   Parameters:
%    - theta is fixed
%    - a Kronecker hypergraph is generated with fixed theta
%    - a random permutation is generated and we attempt to recover the
%    initial permutation
%    - we compare the recovered permutation and initial permutation
%    visually. Note, the 2 aligned graphs do not look identical, but their
%    adjacency matrices are easily isomrphic
%   Functions:
%    - genPermutation.m
%    - permutationProbabilityRatioTest.m
%
% THETA RECOVERY 1 WITH FIXED NODE PERMUTATION
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 19, 2023

%% GRAPH ALIGNMENT 1
clear all; close all; clc;

randSeed = 1; rng(randSeed);                    % Set random seed

kronExp = 3;                                    % Generate kronecker hypergraph
eps = 0.1;
theta = [1 1 1; 0 1 0; 1 0 1];
theta = theta - eps; theta(theta < 0) = eps;
P = theta;
for i=1:kronExp
    P = kron(theta,P);
end
A = (P > 0.25);
[x, y] = find(A == 1); E = [x y]; n = size(A,1);

tic; p = randperm(n);                           % Recover node permutation 
[p, pp] = genPermutation(p, theta, 100000, E);  % using initial theta
t1 = toc;

% Check graph alignment of A into P using the perms
n = size(A,1);
Ap = zeros(n,n);
for i=1:n;    for j=1:n
    Ap(i,j) = A(p(i), p(j));
end; end

figure; 
subplot(1,3,1); imagesc(A); title('Kronecker Graph');
subplot(1,3,2); imagesc(P); title('Kronecker Expansion');
subplot(1,3,3); imagesc(Ap); title('Aligned Graph');

h = figure;
M(size(pp,1)) = struct('cdata',[],'colormap',[]);
for t = 1:size(pp,1)
    Ap = zeros(n,n);
    for i=1:n; for j=1:n
        Ap(i,j) = A(pp(t,i), pp(t,j));
    end; end
    imagesc(Ap);
    M(t) = getframe;
end

%% THETA RECOVERY 1 WITH FIXED NODE PERMUTATION
randSeed = 1; rng(randSeed);                    % Set random seed
clear; clc; close all;                          % Generate random Kronecker graph
theta = [1 1 1; 0 1 0; 1 0 1];
n0 = size(theta, 1);
kExp = 5; n = n0^kExp;
numE = 100 * n;
E = kronGen(theta / sum(theta, 'all'), kExp, numE);
A = zeros(n,n);
for i=1:size(E)
    A(E(i,1),E(i,2)) = 1;
end
[x,y] = find(A == 1); E = [x y];

firstPermItrs = 10000;
itrs = 10000;
theta0 = rand(3,3);
[thetaL, likelihoods] = HyperKronFit("E",E,"theta0",theta0,"maxItrs",20,"firstPermItrs",firstPermItrs,"learningRate",1e-5,"gradSamples",itrs,"v",true,"perm",1:n)

assert(norm(thetaL - theta) < 1e-3)
