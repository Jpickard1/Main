%% Demo
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

% [x, y] = find(A == 1); E = [x y];

clear; clc; close all;
theta = [0.9 0.6
         0.6 0.3];
theta = [1 1 1;
         0 1 0;
         1 0 1];
n0 = size(theta, 1);
kExp = 5; n = n0^kExp;
numE = 100 * n;
E = kronGen(theta / sum(theta, 'all'), kExp, numE)
A = zeros(n,n);
for i=1:size(E)
    A(E(i,1),E(i,2)) = 1;
end
figure; imagesc(A)
[x,y] = find(A == 1); E = [x y];
% filePath = "C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt";
% filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt";
% E = readAdjList(filePath, 0);       % Read file of adjacency list

%% KronFit
theta0 = [0.9 0.6; 0.6 0.3];
eps = 1e-3;
theta0 = theta - eps; theta0(theta0 < 0) = eps;
firstPermItrs = 10000;
itrs = 10000;
theta0 = rand(3,3);
HyperKronFit("E",E,"theta0",theta0,"maxItrs",10,"firstPermItrs",firstPermItrs,"learningRate",1e-5,"gradSamples",itrs,"v",true,"perm",1:n)

%%
n = max(max(E));
A = zeros(n,n);
for i=1:size(E)
    A(E(i,1),E(i,2)) = 1;
end
Ap = zeros(n,n);
for i=1:n
    for j=1:n
        Ap(i,j) = A(p(i), p(j));
    end
end
figure; 
subplot(1,2,1); imagesc(A)
subplot(1,2,2); imagesc(Ap)

%%
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
x = fmincon(@(f, g) sampleGradient(theta0, itrs, firstPermItrs, E), theta0, [], [], [], [], [], [], [], options)

