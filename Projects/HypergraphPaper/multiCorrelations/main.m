%% Multicorrelations
%
%   To argue for studying hyperedge prediction as a function of the method
%   used to construct the observed hypergraph, I construct hypergraph
%   representations of the Mouse Neuron System [Sweeny] using 3 methods of
%   multi-correlation [Drezner, Wang and Zhang, and Taylor], and I compare
%   the structural properties of each.
%
%   REQUIREMENTS:
%   - HAT
%   - Data: Copy_of_mouse 44_fed_fasted_refed traces.csv is a file received
%   from Can Chen that contains time series data for 21 neurons. The file
%   was edited to remove the column headers so that the data are loaded as
%   a matrix.
%
%   OUTPUT:
%       - number of connected components as function of varies threshold
%       - average vertex degree
%       - histograms
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 1, 2023


%% Load Data
clear; close all; clc;
load("Copy_of_mouse 44_fed_fasted_refed traces.csv");
D = Copy_of_mouse_44_fed_fasted_refed_traces;
clear Copy_of_mouse_44_fed_fasted_refed_traces

%% Compute Multi-correlations
k = 5;
[Drezner, idxs] = HAT.multicorrelations(D, k, 'Drezner');
[Wang, ~] = HAT.multicorrelations(D, k, 'Wang');
[Taylor, ~] = HAT.multicorrelations(D, k, 'Taylor');
mkdir(string(k)); % Make output directory

% Most Central Vertices
thresh = 0.9;
T = table();
varNames = ["Drezner", "Wang", "Taylor"];
for i=1:3
    if i==1; M = Drezner; elseif i==2; M = Wang; elseif i==3; M = Taylor; end
    m = find(M > thresh);       % Find idxs in M with value > thresh
    hyperedges = idxs(m, :);    % Sets of vertices with measures > thresh are extracted
    IM = HAT.hyperedges2IM(hyperedges);
    HG = Hypergraph('IM', IM);
    c = HAT.centrality(HG);
    [~, ii] = maxk(c, 5);
    T = [T table(ii,'VariableNames',[varNames(i)])];
end
fName = char(string(k) + "/centrality_" + string(k) + ".txt");
expTable2latex(T, fName);

% Plots

% 3 historgrams of multi-correlations
figure('Position', [10 10 900 250]);
subplot(1,3,1); histogram(Drezner); title('Drezner'); xlabel('Multi-correlation'); ylabel('Number of Hyperedges');
subplot(1,3,2); histogram(Wang); title('Wang and Zheng'); xlabel('Multi-correlation'); ylabel('Number of Hyperedges');
subplot(1,3,3); histogram(Taylor); title('Taylor'); xlabel('Multi-correlation'); ylabel('Number of Hyperedges');
figName = string(k) + "/histograms.png";
saveas(gcf, figName);

% Average vertex degree
figure('Position', [10 10 900 250]);
subplot(1,3,1); title('Degree Distribution');
hold on; xlabel('Multi-correlation Threshold'); ylabel('Average Vertex Degree');
thresholds = 1:-0.05:0;
for i=1:3
    if i==1; M = Drezner; elseif i==2; M = Wang; elseif i==3; M = Taylor; end
    avgVxD = zeros(length(thresholds), 1);
    for ti=1:length(thresholds)
        thresh = thresholds(ti);
        m = find(M > thresh);       % Find idxs in M with value > thresh
        hyperedges = idxs(m, :);    % Sets of vertices with measures > thresh are extracted
        IM = HAT.hyperedges2IM(hyperedges);
        HG = Hypergraph('IM', IM);
        avgVxD(ti) = mean(HG.nodeDegrees);
    end
    plot(thresholds, avgVxD, '-x')
end
legend(["Drezner", "Wang and Zheng", "Taylor"], 'Location','southwest');

% Number of connected components
% figure; 
subplot(1,3,2); title('Hypergraph Connectivity');
hold on; xlabel('Multi-correlation Threshold'); ylabel('Number of Connected Components');
thresholds = 1:-0.005:0.85;
for i=1:3
    if i==1; M = Drezner; elseif i==2; M = Wang; elseif i==3; M = Taylor; end
    numConComp = zeros(length(thresholds), 1);
    for ti=1:length(thresholds)
        thresh = thresholds(ti);
        m = find(M > thresh);       % Find idxs in M with value > thresh
        hyperedges = idxs(m, :);    % Sets of vertices with measures > thresh are extracted
        IM = HAT.hyperedges2IM(hyperedges);
        HG = Hypergraph('IM', IM);
        C = full(HG.cliqueGraph); C = graph(C);
        numConComp(ti) = length(unique(conncomp(C)));
    end
    plot(thresholds, numConComp, '-x')
end

% System entropy
subplot(1,3,3); title('Entropy');
hold on; xlabel('Multi-correlation Threshold'); ylabel('Tensor Entropy');
thresholds = 1:-0.05:0;
for i=1:3
    if i==1; M = Drezner; elseif i==2; M = Wang; elseif i==3; M = Taylor; end
    E = zeros(length(thresholds), 1);
    for ti=1:length(thresholds)
        thresh = thresholds(ti);
        m = find(M > thresh);       % Find idxs in M with value > thresh
        hyperedges = idxs(m, :);    % Sets of vertices with measures > thresh are extracted
        IM = HAT.hyperedges2IM(hyperedges);
        if isempty(IM)
            E(ti) = 0; continue;
        end
        HG = Hypergraph('IM', IM);
        E(ti) = HG.tensorEntropy();
    end
    plot(thresholds, E, '-x')
end
figName = string(k) + "/threeComparisons.png";
saveas(gcf, figName);