% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 15, 2023

% TODO: Set parameters for experiment
numGroups = 10;
% ds = 'Oddysey';
ds = HG;
trials = 1;
predict = true;

O = cell(trials,1); predictions = cell(trials,1); labels=cell(trials,1)
for i=1:trials
    disp('Trial: ' + string(i))
    [O{i}, predictions{i}, labels{i}] = experimentCE(ds, numGroups, predict);
end

% AUC
if predict
    AUC = tableAUC(predictions, labels, numGroups, trials);
    disp(AUC);
end

% Hypergraph Similarity
S = cell(trials,1);
for k=1:trials
    Ok = O{k};
    S{k} = zeros(length(Ok));
    for i=1:length(Ok)
        HGi = Ok{i};
        for j=i+1:length(Ok)
            HGj = Ok{j};
            % S(i,j) = HAT.directSimilarity(HGi,HGj) ;%,"Alpha",) %(HGi.cliqueGraph, HGj.cliqueGraph, 'type', "spanTree");
            S{k}(i,j) = HAT.indirectSimilarity(HGi.cliqueGraph, HGj.cliqueGraph, 'type', "Centrality");
        end
    end
    S{k} = S{k} + S{k}';
end
Smean = zeros(size(S{1}));
Sstd = zeros(size(S{1}));
for i=1:length(Sstd)
    for j=1:length(Sstd)
        s = zeros(length(O), 1);
        for k=1:length(s)
            s(k) = S{k}(i,j);
        end
        Smean(i,j) = mean(s);
        Sstd(i,j) = std(s);
    end
end

[X,Y] = meshgrid(1:numGroups);
%figure; mesh(X,Y,Smean)
%figure; mesh(X,Y,Sstd)
figure; surf(X,Y,Smean)
figure; surf(X,Y,Sstd)

%% Compute statistics on the observed graphs
%{
        true | least popular removed | ... | most popular removed
numV
numE
avdDegree
avgCardinaltiy
entropy
maxA eigenvalue
%}
T = cell(trials,1);
for k=1:trials
    Ok = O{k};
    T{k} = zeros(6, length(Ok));
    for i=1:length(Ok)
        % Get specific hypergraph
        HGki = observedHG(Ok{i});

        % Compute statistics
        [numV, numE]   = size(HGki.IM)
        avgDegree      = mean(sum(HGki,1))
        avgCardinality = mean(sum(HGki,2))
        if strcmp(ds, 'UER')
            entropy    = HAT.tensorEntropy(HGki);
        else
            entropy    = HAT.matrixEntropy(HGki);
        end
        maxEigA        = eig(HGki,1)

        % Save data to table
        T{k}(1,k) = numV;
        T{k}(2,k) = numE;
        T{k}(3,k) = avgDegree;
        T{k}(4,k) = avgCardinality;
        T{k}(5,k) = entropy;
        T{k}(6,k) = magEigA;
    end
    T{k}
end
Tmean = zeros(size(T{1}));
Tstd = zeros(size(T{1}));
for i=1:length(Tstd)
    for j=1:length(Tstd)
        t = zeros(length(O), 1);
        for k=1:length(t)
            t(k) = T{k}(i,j);
        end
        Tmean(i,j) = mean(t);
        Tstd(i,j) = std(t);
    end
end

%%
DS = {'UER', 'DBLP', 'Cora Reference', 'Cora Citation', 'Citeseer Reference', 'Citeseer Citation', 'ArnetMiner Citation', 'Oddysey'};
ds = DS{d};
