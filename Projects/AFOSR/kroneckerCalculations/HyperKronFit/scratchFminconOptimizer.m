%% KronFit with fmincon Optimizer
%
%   This is a test file to determine the feasbility of using fmincon with
%   my current KronFit code.
%
% FUNCTIONS
%   [likelihood, gradient] = sampleGradient(A, theta, itrs, firstPermItrs,
%   debug, directed, E)
%
% VARIABLES
%   theta0: initial conditions
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 16, 2023


theta0 = rand(2,2);
firstPermItrs = 10000;
debug = false;
directed = true;
% [x, y] = find(A == 1); E = [x y];
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% A = full(A(1:2^10,1:2^10));
itrs = 10000;
%%
x = fmincon(@(f, g) sampleGradient(A, theta0, itrs, firstPermItrs, debug, directed, E), theta0, [], [], [], [], [], [], [], options)
x = fmincon(@(f, g) sampleGradient(A, theta, itrs, firstPermItrs, debug, directed, E), theta0, options)
