%% Pairwise to Multiway
%   The purpose of this experiment is to investigate methods of
%   constructing hypergraphs from pairwise data
%
% Auth: Joshua Pickard
% Date: May 12, 2022
%% Set up

close all
clear all
clc

% Generate 1 random network
V = 20;
E = 40;
network = erdos_renyi_network(V, E);
plot(graph(network), 'Layout', 'circle')
title("ER Network (V=" + string(V) +", E=" + string(E) + ")")

