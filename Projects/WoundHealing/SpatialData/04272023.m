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



