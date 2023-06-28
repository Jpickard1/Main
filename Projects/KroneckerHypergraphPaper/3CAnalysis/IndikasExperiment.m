%% Indika's Experiment
%
%   Simulated Polymer Structure --> Synethetic multi-way contacts -->
%   hypergraph representation
%       --> clique expansion
%   kronecker compression (both hypergraphs and clique graphs)
%   analysis
%       --> eigenvalues
%       --> weigheted interactions
%       --> hypergraph similarity
%   
%   What type of claim could I make here:
%       - compressed hypergraphs contain (more/same/less) information than
%       compressed graphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 26, 2023

%% Load multi-way contacts
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 200;
bins = 10;
binSize = 5000 / bins;
BC = cell(nsamps,3);
path2data = "C:\Joshua\MissingData\Projects\KroneckerHypergraphPaper\3CAnalysis\HipHop2HG\";
for s=1:3
    for i=1:nsamps
        filename = "adjList_" + SIMPARMS(s) + "_" + string(i) + "_0.01.txt";
        filepath = path2data + filename; disp(filepath);
        E = readmatrix(filepath);
        B = zeros(bins, bins, bins);
        for e=1:size(E,1)
            edge = E(e,:);
            r = range(ceil(edge/binSize)*binSize);
            % Internal edge
            if r >= binSize
                % disp(edge); 
                % disp(ceil(edge/binSize)*binSize);
                % disp(r);
                edge = floor((edge - 1) / binSize) + 1;
                pedge = perms(edge);
                pedge = unique(pedge,'rows');
                for p=1:size(pedge,1)
                    B(pedge(p,1), pedge(p,2), pedge(p,3)) = B(pedge(p,1), pedge(p,2), pedge(p,3)) + 1;
                end
            end
        end
        BC{i,s} = B;
    end
end

% Perform Clique Expansion
BCC = cell(nsamps,3);
for s=1:3
    for i=1:nsamps
        B = BC{i,s};
        C = tensorprod(B, ones(bins,1), 1, 1);
        BCC{i,s} = C;
    end
end

% save("BC_nsamps_200_binSize_5.mat", "BC");
% save("BC_nsamps_200_binSize_10.mat", "BC");
% clear
% load("BC_nsamps_200_binSize_10.mat"); bins = 10; nsamps = 10;
% SIMPARMS = ["OFF","ON","HIGH"];
%% Construct Laplacians
LC = cell(size(BC));
LCC = cell(size(BCC));
for s=1:3
    for i=1:nsamps
        A = BC{i,s};
        D = zeros(5,5,5);
        AS1 = sum(A,1); AS2 = sum(AS1, 2); Asum = reshape(AS2, [5 1]);
        for j=1:5
            D(j,j,j) = Asum(j);
        end
        L = D - A;
        LC{i,s} = L;
        A = BCC{i,s};
        L = diag(sum(A)) - A;
        LCC{i,s} = L;
    end
end

%% Degree sequence
DSC = cell(3,1);
DSCC = cell(3,1);
figure;
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    degreeSequence = zeros(nsamps, 5);
    if si <= 3
        for i=1:nsamps
            B = BC{i,s};
            degreeSequence(i,:) = reshape(sum(sum(B)), [1 5]);
        end
        DSC{s} = degreeSequence;
    else
        for i=1:nsamps
            B = BCC{i,s};
            degreeSequence(i,:) = reshape(sum(B), [1 5]);
        end
        DSCC{s} = degreeSequence;
    end
    subplot(2,3,si);
    boxplot(degreeSequence);
    title(SIMPARMS(s));
    if s == 1
        ylabel('Degree');
    end
    xlabel('Polymer Compartment');
    disp(mean(degreeSequence));
end

%% calculate eigenvalues of A
eigA = cell(3,1);
eigAC = cell(3,1);
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    degreeSequence = zeros(nsamps, 5);
    if si <= 3
        for i=1:nsamps
            B = BC{i,s};
            eigA{i,s} = heig(B);
        end
    else
        for i=1:nsamps
            B = BCC{i,s};
            eigAC{i,s} = eig(B);
        end
    end
end

%% calculate eigenvalues of L
eigL = cell(3,1);
eigLC = cell(3,1);
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    degreeSequence = zeros(nsamps, 5);
    if si <= 3
        for i=1:nsamps
            L = LC{i,s};
            eigL{i,s} = heig(L);
        end
    else
        for i=1:nsamps
            L = LCC{i,s};
            eigLC{i,s} = eig(L);
        end
    end
end


%% Plot eigenvalue distributions
figure;
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    subplot(2,3,si);
    if si <= 3
        eigValues = [];
        for i=1:nsamps
            eigValues = [eigValues, eigA{i,s}];
        end
    else
        eigValues = [];
        for i=1:nsamps
            eigValues = [eigValues, eigAC{i,s}];
        end
    end
    histogram(eigValues);
    title(SIMPARMS(s));
    if s == 1
        ylabel('Count');
    end
    xlabel('Eigenvalues');
end

figure;
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    subplot(2,3,si);
    if si <= 3
        eigValues = zeros(nsamps, 1);
        for i=1:nsamps
            eigValues(i) = max(eigA{i,s});
        end
    else
        eigValues = zeros(nsamps, 1);
        for i=1:nsamps
            eigValues(i) = max(eigAC{i,s});
        end
    end
    boxplot(eigValues);
    % historgram(eigValues);
    title(SIMPARMS(s));
    if s == 1
        ylabel('Count');
    end
    xlabel('Eigenvalues');
end

figure;
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    subplot(2,3,si);
    if si <= 3
        eigValues = [];
        for i=1:nsamps
            eigValues = [eigValues, eigL{i,s}];
        end
    else
        eigValues = [];
        for i=1:nsamps
            eigValues = [eigValues, eigLC{i,s}];
        end
    end
    histogram(eigValues);
    title(SIMPARMS(s));
    if s == 1
        ylabel('Count');
    end
    xlabel('Eigenvalues');
end

figure;
for si=1:6
    s = mod(si-1, 3)+1; disp(si); disp(s)
    subplot(2,3,si);
    if si <= 3
        eigValues = zeros(nsamps, 1);
        for i=1:nsamps
            eigValues(i) = max(mink(eigA{i,s}, 2));
        end
    else
        eigValues = zeros(nsamps, 1);
        for i=1:nsamps
            disp(eigAC{i,s});
            ll = mink(eigAC{i,s},2);
            disp(ll');
            eigValues(i) = max(ll);
        end
    end
    boxplot(eigValues);
    % historgram(eigValues);
    title(SIMPARMS(s));
    if s == 1
        ylabel('Count');
    end
    xlabel('Eigenvalues');
end

%% Hypergraph construction
%   Threshold to consider only the top 5 most significant hyperedges
HGC = cell(size(BC));
hedges = nchoosek(1:bins,3);
for s=1:3
    disp(s);
    for i=1:nsamps
        B = BC{i,s};
        hedgeWeights = zeros(size(hedges, 1), 1);
        for j=1:size(hedges,1)
            hedgeWeights(j) = B(hedges(j,1),hedges(j,2), hedges(j,3));
        end
        % [~, idxs] = find(hedgeWeights ./ max(hedgeWeights) > 0.5);
        [~, idxs] = maxk(hedgeWeights, round(0.5*size(hedges,1)));
        HGedges = hedges(idxs,:);
        HGC{i,s} = Hypergraph('IM', HAT.hyperedges2IM(HGedges));
        disp("===========")
        disp(HGedges);
    end
end

GC = cell(size(BC));
edges = nchoosek(1:bins,2);
for s=1:3
    disp(s);
    for i=1:nsamps
        B = BCC{i,s};
        edgeWeights = zeros(size(edges, 1), 1);
        for j=1:size(edges,1)
            edgeWeights(j) = B(edges(j,1),edges(j,2));
        end
        % [~, idxs] = find(edgeWeights ./ max(edgeWeights) > 0.5);
        [~, idxs] = maxk(edgeWeights, round(0.5*size(edges,1)));
        Gedges = edges(idxs,:);
        GC{i,s} = Hypergraph('IM', HAT.hyperedges2IM(Gedges));
        disp("===========")
        disp(Gedges);
    end
end

%% Hypergraph entropy

HGE = zeros(size(BC));
GE = zeros(size(BC));
for s=1:3
    for i=1:nsamps
        HGE(i,s) = HAT.tensorEntropy(HGC{i,s});
        GE(i,s) = HAT.tensorEntropy(GC{i,s});
    end
end

figure;
subplot(1,2,1); boxplot(HGE); xticklabels(SIMPARMS); title('Hypergraph Entropy');
subplot(1,2,2); boxplot(GE); xticklabels(SIMPARMS); title('Graph Entropy');

figure;
for si=1:6
    s = mod(si-1, 3) + 1; disp(si); disp(s);
    subplot(2,3,si);
    if si <= 3
        etp = HGE(:,s);
    else
        etp = GE(:,s);
    end
    histogram(etp,10);
end

%% Hypergraph Similarity
%   Compute the distance matrix between all hypergraphs

% Create all adjacency tensors
AC = cell(size(HGC));
for s=1:3
    for i=1:nsamps
        HG = HGC{i,s};
        AC{i,s} = HG.adjTensor;
    end
end

%% Hypergraph hamming similarity
D = zeros(3*nsamps, 3*nsamps);
for idx1=1:(3*nsamps)
    s1 = ceil(idx1 / nsamps);
    i1 = idx1 - ((s1-1) * nsamps);
    A1 = AC{i1, s1};
    disp(idx1);
    disp(string(s1) + ", " + string(i1));
    for idx2=idx1:(3*nsamps)
        s2 = ceil(idx2 / nsamps);
        i2 = idx2 - ((s2-1) * nsamps);
        A2 = AC{i2, s2};
        d = sum(abs(A1-A2), 'all')/(10^3-10);
        D(idx1, idx2) = d;
    end
end

D = D + D';
figure; imagesc(D)

Z = linkage(D);
figure;
dendrogram(Z)



