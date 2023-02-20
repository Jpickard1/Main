%% Preamble

clear all;
clc;
% In DGC directory
load('data.mat')
% Select the important columns
% Columns:
% - 322:328 Stages
% - 329     Shae-Kelly Human Islet
% - 320, 321, 319   Other islet
D = data(:,[322:328 329 320 321 319]);
D = D{:,:};
D(isnan(D) | isinf(D)) = 0;

% genes where 7 or more are 0
%   Reduces 22083 to 18706
for i=1:length(D)
    if length(find(D(i,:) == 0)) > 6
        D(i,:) = [];
    end
end

numCells = min(size(D));
for i=1:numCells
    D(:,i) = D(:,i) / sum(D(:,i));
end
sum(D)

%% Covariance matrix
C = corrcov(cov(D));
figure; h = HeatMap(C, 'XData', ["S0","S1","S2","S3","S4","S5","S6","Islet (Shae)","Islet (Jens)","Islet (Michael)","Islet (Anthony)"], 'YData', ...
        ["S0","S1","S2","S3","S4","S5","S6","Islet (Shae)","Islet (Jens)","Islet (Michael)","Islet (Anthony)"])
h.Annotate = true;
title('Correlation Matrix');
saveas(gcf, 'correlation.png')

C = cov(D);
figure; heatmap(C, 'XData', ["S0","S1","S2","S3","S4","S5","S6","Islet (Shae)","Islet (Jens)","Islet (Michael)","Islet (Anthony)"], 'YData', ...
        ["S0","S1","S2","S3","S4","S5","S6","Islet (Shae)","Islet (Jens)","Islet (Michael)","Islet (Anthony)"])
title('Covariance Matrix');
saveas(gcf, 'Covariance.png')



%% Make Distance matrix

distances = ["euclidean", "squaredeuclidean", "cityblock", "minkowski", "chebychev", "mahalanobis", "cosine", "correlation", "spearman", "hamming", "jaccard"];
for i=1:length(distances)
    d = squareform(pdist(D', distances(i)));
    figure;
    heatmap(d, 'XData', ["S0","S1","S2","S3","S4","S5","S6","Islet (Shae)","Islet (Jens)","Islet (Michael)","Islet (Anthony)"], 'YData', ...
        ["S0","S1","S2","S3","S4","S5","S6","Islet (Shae)","Islet (Jens)","Islet (Michael)","Islet (Anthony)"])
    title(distances(i) + ' Distance')
    saveas(gcf, distances(i) + '_distance.png')
end

%{
d = squareform(pdist(D', 'cityblock'));
figure;
heatmap(d, 'XData', ["Col 1","Col 2","Col 3","Col 4","Col 5","Col 6","Col 7","Col 8","Islet Cells","Fibroblast"], 'YData', ...
    ["Col 1","Col 2","Col 3","Col 4","Col 5","Col 6","Col 7","Col 9","Islet Cells","Fibroblast"])
title('Cityblock Distance')
saveas(gcf, 'cityblock_distance.png')


% ax = gca;
% ax.YData = {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}, {'g'}, {'h'}, {'i'}, {'j'}}
xticks(10:19)
length(find(isnan(D)))

length(find(isnan(D)))

D = data(:,[319:326, 328]);
D = D{:,:};
for i=1:9
    D(:,i) = D(:,i) / sum(D(:,i));
end
sum(D)

t = array2table(D, 'VariableNames',{'KS 1 - s1', 'KS 2 - s1', 'KS 3 - s1', 'KS 4 - s1', 'KS 5 - s1', 'KS 6 - s1', 'KS 7 - s1', 'KS 8 - s1', 'Islet s1'});
data = [data t];
save data.mat data
%}

%%
% Remove old data from them
% data(:,329:337) = [];
% data(:,319:326) = [];

