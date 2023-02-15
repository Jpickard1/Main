% DRIVERLP (Cold Ends)
%
%   This driver follows the following workflow:
%       1. Set experiment parameters
%       2. Load hypergraph
%       3. Create observed systems according to exact cold ends problem
%       4. Hyperedge Prediction
%       5. View AUC of each prediction method on each missing probe set
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023

% TODO: Set parameters for experiment
numGroups = 4;
% thresh = 0.20;
% itrs = 1;
% trials = 1;

% TODO: Set dataset
d = 3;
DS = {'DBLP', 'Cora Reference', 'Cora Citation', 'Citeseer Reference', 'Citeseer Citation', 'ArnetMiner Citation', 'Oddysey'};
ds = DS{d};

% TODO: Set observation bias

% Load data as hypergraph
HGt = load_HG(ds);

% Observations
% Generate observations of hypergraph system
[O, E] = observeSystemCE(HGt,numGroups);
for i=1:length(O)
    disp(size(observedHG(O{i}).IM))
end

% figure; plot(1:length(E),E)
% Prediction
% Make false hyperedges
F = falseHyperedges(HGt, 2*size(HGt.IM,2));

% Hyperedge prediction on each observed hypergraph
predictions = cell(numGroups,1);
labels = cell(numGroups,1);
for i=1:numGroups
    HGobserved = O{i};
    disp(size(observedHG(HGobserved).IM));
    [predictions{i}, labels{i}] = predictHyperedges(HGt, HGobserved, F);
end

% AUC and F1 Scores
t = [];
for i=1:length(O) % (num observations)
    c = zeros(7,2);
    for j=1:7 % (num predictors)
        auc = zeros(1,trials);
        for k=1:trials
            p = predictions{i}(:,j);
            l = labels{i};
            [~, ~, ~, auc(k)] = perfcurve(l, p, 1);
        end
        c(j,1) = mean(auc);
        c(j,2) = std(auc);
    end
    t = [t c];
end
disp(t);

