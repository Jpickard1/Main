function [T] = scratchPrelimExperiment2(pKnown)
%DEMOEXPERIMENT This is a scratch file to implement my experiment from
% March 10, 2023 quickly for multiple iterations. See 
% BLinkPrediction/scratch.m for a summary of whats happening here.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 13, 2023

F = {@adamicAdar_index, @commonNeighbors_index, ...
    @jaccard_index, @leichtHolmeNewman_index }; %, ...
    % @averageCommuteTime_index, @randomWalkWithRestart_index};

% Load data
load('C:\Users\picka\Documents\my_projects\DBTM\Main\Data\USAir\adjacency_matrix.mat')
A = A + A'; A = (A>0); A = real(A); n = size(A,1);
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
B = [400 500 size(Eu, 1)];

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

end

