% April 27, 2023
%
%   This file constructs hypergraphs from Cooper's spatial data.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 27, 2023

clear; close all; clc;

addpath(genpath("/mnt/c/Joshua/Main/Projects/AFOSR/"));
addpath(genpath("/mnt/c/Joshua/Hypergraph-Analysis-Toolbox"));
load("C:\Joshua\MissingData\Projects\WoundHealing\SpatialData\04272023_Spatial_HGs_1.mat");

AC = C.adjTensor;
AG = G.adjTensor;

% Compute Kronecker Hypergraph
K = superkron(AC, AG);
V = kronString(CV, GV);

% Get the linear indices of all nonzero elements
idx = find(K);
ES = [];
[ES(:,1), ES(:,2), ES(:,3)] = ind2sub(size(K), idx);
IM = HAT.hyperedge2IM(ES);
IM = unique(IM', 'rows')';

ES = zeros(size(IM,2), 3);
for i=1:size(IM,2)
    ES(i,:) = find(IM(:,i));
end

csvwrite('04272023_edgeset.csv', ES);

% Open file for writing
fid = fopen('04272023_vertexLabels.txt', 'w');
% Write each string in the array to a separate line in the file
for i = 1:numel(S)
    fprintf(fid, '%s\n', S(i));
end
% Close file
fclose(fid);
