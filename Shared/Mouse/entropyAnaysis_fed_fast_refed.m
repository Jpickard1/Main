% Read the data
T = readtable('mouse 44_fed_fasted_refed traces.csv');
data = table2array(T(1:end, 2:end)); 
time = table2array(T(1:end, 1));

n = size(data, 2); % number of neurons
stageStart = [0, 8.57002e+04, 1.7769786e+05]; % starting point for each stage
stageStartIdx = find(ismember(time, stageStart)>0);

% Find the time series data for each stage
fedStage = data(stageStartIdx(1):stageStartIdx(2)-1, :);
fastStage = data(stageStartIdx(2):stageStartIdx(3)-1, :);
refedStage = data(stageStartIdx(3):end, :);

% Correlation analysis for each stage
stageMatrix = {fedStage, fastStage, refedStage};
entropy = zeros(3, 1);
rawCorrMatrix = zeros(n, n, 3);
degree_centrality = zeros(n, 3);
eps = 0.7; % epsilon
for i = 1:3
    rawCorrMatrix(:, :, i) = corrcoef(stageMatrix{i});
    corrMatrix = abs(corrcoef(stageMatrix{i}) - eye(n)); 
    corrMatrix(corrMatrix>=eps) = 1; % binarize
    corrMatrix(corrMatrix<eps) = 0; % binarize
    
    % Construct the graph
    G = graph(corrMatrix);
    L = laplacian(G); 
    
    % Compute the eigenvalues
    eigVal = eig(full(L));
    nonzeroEigVal = eigVal(eigVal>1e-8);
    normalizedEig = nonzeroEigVal/sum(nonzeroEigVal);
    entropy(i) = -sum(normalizedEig.*log(normalizedEig));
    
    % Compute degree centrality
    degree_centrality(:, i) = centrality(G,'degree');
end


% Correlation matrix
figure, 
for i = 1:3
    subplot(1, 3, i)
    imagesc(rawCorrMatrix(:, : ,i));
    axis square
    colormap hot
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
end

% Bar plot the entropy and centrality
figure, 
subplot(1, 2, 1)
bar(entropy)
axis square
ylim([1.5 2.5])
subplot(1, 2, 2)
bar(degree_centrality, 'stack')
axis square
xlim([0 21])