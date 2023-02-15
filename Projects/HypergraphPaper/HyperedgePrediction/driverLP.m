% DRIVERLP (Link Prediction)
%
%   This driver follows the following workflow:
%       1. Set experiment parameters
%       2. Load hypergraph
%       3. Observe hypergraphs with varying biases
%       4. Hyperedge Prediction
%       5. View AUC of each prediction method on each type of bias
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

% TODO: Set parameters for experiment
thresh = 0.20;
itrs = 1;
trials = 1;

% TODO: Set dataset
d = 4;
DS = {'DBLP', 'Cora Reference', 'Cora Citation', 'Citeseer Reference', 'Citeseer Citation', 'ArnetMiner Citation', 'Oddysey'};
ds = DS{d};

% TODO: Set observation bias
rL = -0.5; rH = 0.5;
r = [rL, 0, rH];

% Load data as hypergraph
HGt = load_HG(ds);

% Observations
% Generate observations of hypergraph system
O = cell(trials,1);
for t=1:trials
    O{t} = observeSystemLP(HGt, thresh, r);
end

for t=1:trials
    for i=1:length(O{t})
        disp(size(observedHG(O{t}{i}).IM))
    end
end

% Prediction

% Hyperedge prediction on each observed hypergraph
predictions = cell(length(O), 6);
labels = cell(length(O), 6);
for t=1:trials
    for i=1:6
        HGobserved = O{t}{i};
        disp(size(observedHG(HGobserved).IM));
        [predictions{t, i}, labels{t,i}] = predictHyperedges(HGt, HGobserved);
    end
end

%% AUC and F1 Scores
t = [];
for i=1:6
    c = zeros(7,2);
    for j=1:7
        auc = zeros(1,trials);
        for k=1:trials
            p = predictions{k,i}(:,j);
            l = labels{k,i};
            [~, ~, ~, auc(k)] = perfcurve(l, p, 1);
        end
        c(j,1) = mean(auc);
        c(j,2) = std(auc);
    end
    t = [t c];
end
disp(t);
% tt = t;
%% Plots
% Score hyperedges
HOMs = ["OVH","OVN","OVL","OEH","OEN","OEL"];
vxSim = ["CN", "LHN", "SA", "SO", "HP", "HD", "JC"];
nbars = 10;
figure;
for i=1:6
    for j=1:7
        subplot(6,7,((i-1)*7) + j);
        p = predictions{i}(:,j);
        d1 = p(1:floor(end/2));
        d2 = p(floor(end/2)+1:end);
        histogram(d1, nbars); hold on; histogram(d2, nbars);        
        if j == 1
            ylabel(HOMs(i));
        end
        if i==6
            xlabel(vxSim(j));
        end
    end
end

HOMs = ["OVH","OVN","OVL","OEH","OEN","OEL"];
vxSim = ["CN", "LHN", "SA", "SO", "HP", "HD", "JC"];
figure;
for i=1:6
    for j=1:7
        subplot(6,7,((i-1)*7) + j);
        p = predictions{i}(:,j);
        [x, y, ~, auc] = perfcurve(labels{i}, p, 1);
        plot(x, y); title(string(auc));
        if j == 1
            ylabel(HOMs(i));
        end
        if i==6
            xlabel(vxSim(j));
        end
    end
end

% figure; perfcurve(labels{1}, p, 1)

%% figure; histogram(rand(1,1000), 20); hold on; histogram(rand(1,1000), 20);


