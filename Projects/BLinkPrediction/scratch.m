%% Batched Prediction
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 9, 2023

%% Marck 10, 2023


%% Demo Experiment
%   Here I partition the edge set into observed and (n) unobserved edges. I
%   generate n false or negative set of edges. Then for a few link
%   prediction functions F, I perdict a total of n edges. The n edges are
%   predicted in groups of size b, and after being predicted the edges are
%   assumed true when predicting the next batch of edges. I vary the value
%   of b, and we see that the performance of each prediction function
%   varies significantly with respect to the batch size b.
%
%   POSSIBLE OUTPUT from ">> disp(T)"
%      0.76176     0.77273     0.80564    0.82602     0.8511    0.88715    0.84483    0.93103     0.76019
%      0.77743      0.7884     0.82445    0.85737    0.86364    0.90125    0.87774    0.92633     0.78997
%    0.0031348    0.018809    0.097179    0.17555    0.25392    0.41066    0.25392     0.5674    0.012539
%    0.0031348    0.018809    0.097179    0.17555    0.25392    0.41066    0.25392     0.5674    0.012539
%    0.0031348    0.018809    0.097179    0.17555    0.25392    0.41066    0.25392     0.5674    0.012539
%      0.18025     0.21003     0.25862    0.30878     0.3558    0.47335     0.3558    0.60345     0.19749
%    0.0031348    0.018809    0.097179    0.17555    0.25392    0.41066    0.25392     0.5674    0.012539
%    0.0031348    0.018809    0.097179    0.17555    0.25392    0.41066    0.25392     0.5674    0.012539
%
clear; close all; clc

% TODO: Set parameters
pKnown = 0.7;   % number of edges to observe

% Set list of link prediction functions
F = {@adamicAdar_index, @commonNeighbors_index, @hubDepressed_index, @hubPromoted_index, ...
    @jaccard_index, @leichtHolmeNewman_index, @salton_index, @sorensen_index }; %, ...
    % @averageCommuteTime_index, @randomWalkWithRestart_index};

% Load data
load('C:\Users\picka\Documents\my_projects\DBTM\Main\Data\USAir\adjacency_matrix.mat')
A = (A>0); A = real(A); n = size(A,1);
A = A - tril(A);
[vi, vj] = find(A == 1);
E = [vi vj];                    % Get edge set

% 2. Remove edges
E = E(randperm(size(E, 1)), :); % Randomly permute edge set

Ek = E(1:round(size(E,1) * pKnown), :);
Eu = E(round(size(E,1) * pKnown) + 1:end, :);

% Get edges to predict
Ec = nchoosek(1:n, 2);         % list all possible edges
Ec = setdiff(Ec, E, 'rows');   % remove edges that already exist
Ec = Ec(randperm(size(Ec, 1)), :);
Ep = [Ec(1:size(Eu), :); Eu];

% Set batch sizes
B = [20 50 100 150 200 300 400 500 size(Eu, 1)];

% round(exp(1:1:log(size(Eu, 1))))
% B(end + 1) = size(Eu, 1);

T = table();
for bi=1:length(B)
    b = B(bi);
    disp(b);
    bT = cell(length(F), 1);
    for fi=1:length(F)
        f = F{fi};
        disp(f);
        Ei = bpredict(n, Ek, size(Eu, 1), b, f, Ep);
        acc = (size(Eu, 1) - size(setdiff(Eu, Ei, 'rows'), 1)) / size(Eu, 1); % / size(Eu, 1);
        T{fi, bi} = acc;
    end
end
disp(T)

%% Chicago (118) and DTW (112) neighborhoods
clear
load('C:\Users\picka\Documents\my_projects\DBTM\Main\Data\USAir\adjacency_matrix.mat')
A = (A>0); A = real(A); n = size(A,1);
nc = find(A(118,:) == 1);
nd = find(A(112,:) == 1);

Ad = A(nd, nd);
figure; plot(graph(Ad))


%% March 9, 2023

clear; close all; clc

% TODO: Parameters
n = 100;        % network size
pKnown = 0.9;   % number of edges to observe
b = 2;          % batch size

% 1. Make network
A = erdos_renyi_network(n, round(0.5*nchoosek(n,2))); A = real(A);
A = A - tril(A);
[vi, vj] = find(A == 1);
E = [vi vj];                    % Get edge set

% 2. Remove edges
E = E(randperm(size(E, 1)), :); % Randomly permute edge set

Ek = E(1:round(size(E,1) * pKnown), :);
Eu = E(round(size(E,1) * pKnown) + 1:end, :);
% size(Eu,1) + size(Ek,1)
% eig

% 3. call link prediction function
Ei = bpredict(n, Ek, size(Eu, 1), b, @commonNeighbors_index);

disp(size(setdiff(Eu, Ei, 'rows'), 1) / size(Eu, 1))
