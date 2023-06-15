%% Test Analysis
%
%   This file checks the output of test data from Great Lakes.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% June 14, 2023

clear; clc; close all;

load("syntheticData1.mat")
figure; plot(likelihoods)

load("syntheticData2.mat")
figure; plot(likelihoods)
