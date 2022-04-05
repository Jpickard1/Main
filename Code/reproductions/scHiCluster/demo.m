%% Load
% readtable('chr19.csv');
csvread('chr19.csv');
A = (ans(2:end,2:end));
% A.data
%% 
mat = driver_scHiCluster(A)