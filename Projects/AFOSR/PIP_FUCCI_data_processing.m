%% PIP FUCCI DATA PROCESSING
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 2, 2023

clear; clc; close all;

T = readtable('keggCellCycleGenes.csv');
geneNames = T.GeneName;

% Load DGC data to find multicorrelations in these genes
load('C:\Users\picka\Documents\my_projects\DBTM\DataGuidedControl\Data\data.mat')

indexGeneNames = data.Properties.RowNames;

[~, ia, ~] = intersect(indexGeneNames, geneNames);

data = removevars(data,{'external_gene_name'});
D = data{ia,:}';

D(isnan(D)) = 0;
[M, idxs] = HAT.multicorrelations(D, 3, 'Drezner');

figure; histogram(M); xlabel('Multi-correlation'); ylabel('Number of Hyperedges'); title('KEGG Cell Cycle Multi-Correlations')

%%

thresh = 0.998;
m = find(M > thresh);       % Find idxs in M with value > thresh
hyperedges = idxs(m, :);     % Sets of vertices with measures > thresh are extracted
IM = HAT.hyperedges2IM(hyperedges);
HG = Hypergraph('IM', IM);

C = graph(full(HG.cliqueGraph));
disp(length(unique(conncomp(C))));

figure; plot(C, 'Layout', 'circle'); title('Clique Expanded Cell Cycle Graph');
figure; HG.plot(); title('KEGG Cell Cycle Hypergraph'); xlabel('Hyperedges'); ylabel('Vertices');

