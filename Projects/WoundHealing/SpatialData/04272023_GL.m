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


% April 27, 2023
%
%   This file constructs hypergraphs from Cooper's spatial data.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 27, 2023

%% Great Lakes Preamble
addpath(genpath("/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/"));
% addpath(genpath("/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/"));
addpath(genpath("/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox"));

%% Load Data
dataPath = "/nfs/turbo/umms-indikar/shared/projects/spatial_transcriptomics/data/kronecker_stuff/";
Atab = readtable(dataPath + "A.csv");
Btab = readtable(dataPath + "B.csv");
coordstabl = readtable(dataPath + "coords.csv");

%% Cell Type Hypergraph
%   Vertex set is the cell types (n=11). Edge sets are groups of 3 cells
%   that have preference for interacting with one another.

% Cooper's data --> incidence matrix
IMw = Atab{:,2:end};            % Extract numerical data
[~, ES] = maxk(IMw, 3, 1);      % Find 3 most predicted cells at each spot
IM = HAT.hyperedges2IM(ES');    % Convert edge set to incidence matrix
IM = unique(IM', 'rows')';      % Remove duplicate hyperedges

% Set Hypergraph
C = Hypergraph('IM', IM);       % Create hypergraph object
CV = Atab{:,1};                 % Save vertex set of hypergraph

%% Transcript Hypergraph
IMw = Btab{:,2:end};            % Extract numerical data
[~, ES] = maxk(IMw, 3, 1);      % Find 3 most predicted genes at each spot
IM = HAT.hyperedges2IM(ES');    % Convert edge set to incidence matrix
IM = unique(IM', 'rows')';       % Remove duplicate hyperedges

vxDegree = sum(IM,1);           % Get degree of vertices
[~, vxKeep] = find(vxDegree~=0);% Get vertices with nonzero degree
IM = IM(vxKeep, :);           % Reduce the vertex set

% Set Hypergraph 
G = Hypergraph('IM', IM);       % Create hypergraph object
GV = Btab{vxKeep,1};            % Save vertex set of hypergraph

%% Save hypergraphs
save('/nfs/turbo/umms-indikar/Joshua/Main/Projects/WoundHealing/Sparial04272023_Spatial_HGs_1.mat', 'G','C','GV','CV');
