% k Uniform Neurons 3 Ways
%
%   Here I take the mouse feeding neuron data and construct a hypergraph
%   using all 3 multi-correlations currently implemented in HAT. Then, I
%   compate the properties of each hypergraph.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 15, 2023

%% Load data and build hypergraphs
%   the threshold for all hypergraphs is set to be the minimal threshold
%   that keeps the network connected

% Load data
D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');

% Drezner
[M, idxs] = HAT.multicorrelations(D, 3, 'Drezner');
% figure; histogram(M)
t = 0.83; % Disconnected at 0.84
hyperedges = idxs(find(M>t),:);
IM = HAT.hyperedge2IM(hyperedges);
HG = Hypergraph('IM', IM);
% figure; plot(graph(HG.cliqueGraph)) % To check the graph is fully connected
HGd = HG;

% Wang and Zheng
[M, idxs] = HAT.multicorrelations(D, 3, 'Wang');
% figure; histogram(M)
t = 0.83; % Disconnected at 0.84
hyperedges = idxs(find(M>t),:);
IM = HAT.hyperedge2IM(hyperedges);
HG = Hypergraph('IM', IM);
% figure; plot(graph(HG.cliqueGraph)) % To check the graph is fully connected
HGw = HG;

% Taylor
[M, idxs] = HAT.multicorrelations(D, 3, 'Taylor');
% figure; histogram(M)
t = 0.92; % Disconnected at 0.93
hyperedges = idxs(find(M>t),:);
IM = HAT.hyperedge2IM(hyperedges);
HG = Hypergraph('IM', IM);
% figure; plot(graph(HG.cliqueGraph)) % To check the graph is fully connected
HGt = HG;

HG = {HGd, HGt, HGw};
clearvars -except HG

%% Compate hypergraph structure

S = zeros(3,3);
for i=1:3
    for j=i+1:3
        S(i,j) = HAT.directSimilarity(HG{i},HG{j});
    end
end
S = S + S';
