function [O, predictions, labels] = experimentCE(ds, numGroups, predict)
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

if nargin == 2
    predict = true;
end

% Load data as hypergraph
if strcmp(class(ds), 'Hypergraph')
    HGt = ds;
else
    HGt = load_HG(ds);
end

% Make false hyperedges
F = falseHyperedges(HGt, size(HGt.IM,2));

% Observations
% Generate observations of hypergraph system
[O, E] = observeSystemCE(HGt,numGroups);
for i=1:length(O)
    disp(size(observedHG(O{i}).IM))
end

if ~predict
    predictions = cell(numGroups,1);
    labels = cell(numGroups,1);
    return
end

% figure; plot(1:length(E),E)
% Prediction
% Hyperedge prediction on each observed hypergraph
for i=1:numGroups
    HGobserved = O{i};
    % disp(size(observedHG(HGobserved).IM));
    [predictions{i}, labels{i}] = predictHyperedges(HGt, HGobserved, F);
end

end


