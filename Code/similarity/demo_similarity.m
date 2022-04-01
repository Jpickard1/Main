% Description: This file is a driver file that computes missing links based
% using different similarity indices.

%% Load data
close all;
clear all;
clc;
data_path = data_path();
data_set = "USAir";
load(data_path + data_set + "/adjacency_matrix.mat");

%% Split data
[known_network, unknown_network] = edge_removal(A, 'R');

%% Compute similarity matrix
similarity = compute_similarity(known_network, 'AA');
predicted_network = similarity_predict(known_network, similarity, 0.5);

max(max(predicted_network))
sum(sum(predicted_network))
sum(sum(unknown_network))

%% Measure performance

AUC = eval_metrics(known_network, unknown_network, predicted_network);
disp(AUC)
