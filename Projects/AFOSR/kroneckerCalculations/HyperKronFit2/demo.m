%% Demo
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

% [x, y] = find(A == 1); E = [x y];

clear; clc; close all;
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt";
E = readAdjList(filePath, 0);       % Read file of adjacency list
theta0 = rand(2,2);
firstPermItrs = 10000;
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
itrs = 10000;

x = fmincon(@(f, g) sampleGradient(theta0, itrs, firstPermItrs, E), theta0, [], [], [], [], [], [], [], options)

