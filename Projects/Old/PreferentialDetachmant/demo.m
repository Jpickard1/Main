%% Preferential Detachment Algorithm (PDA)
%   The Preferential Detachment Algorithm (PDA) is motivated by the 
%   Barabasi-Albert model of Preferential Attachment scale free networks. 
%   The PA algorithm results in scale free networks where the degree
%   distribution in the network obeys a power law. This process works by
%   adding edges with preference towards the most central nodes (degree 
%   centrality although other measures work as well).
%
%   I would like to engineer a reversed process to deconstruct a graph into
%   multiple subgraphs that have similar characteristics of the original
%   graph.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2022

%% Preamble
close all
clear
%%

% A = generate_random_network('ER', 20, 0.5);
A = generate_random_network('SF', 50, 0.5);
issymmetric(A)
cords = circleCords(A);
figure; plot(graph(A), 'XData', cords(:,1), 'YData', cords(:,2))

% figure; plot(sort(sum(A)))
% figure; gplot(A, cords)

%% Scratch

A = [0 1 1 1;
     1 0 1 0;
     1 1 0 0;]

%% 
n = 1000;
A = true(n,n);
for i=1:n; A(i,i) = false; end

while sum(sum(A)) >= n-1
    degree = sum(A);
    
    degree = degree ./ sum(degree);
    
end





