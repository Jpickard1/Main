function [t] = tableAUC(predictions, labels, numCols, trials)
%TABLEAUC Summary of this function goes here
%   Detailed explanation goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 15, 2023

%%
t = [];
for i=1:numCols % (num observations)
    c = zeros(7,2);
    for j=1:7 % (num predictors)
        auc = zeros(1,trials);
        for k=1:trials
            p = predictions{k}{i}(:,j);
            l = labels{k}{i};
            [~, ~, ~, auc(k)] = perfcurve(l, p, 1);
        end
        c(j,1) = mean(auc);
        c(j,2) = std(auc);
    end
    t = [t c];
end
end

