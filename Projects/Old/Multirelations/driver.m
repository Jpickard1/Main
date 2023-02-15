%% Multirelations
%
% PROCESS:
%   1. construct random hypergraph incidence matrix
%   2. clique expand to generate graph data
%   3. reconstruct hypergraph using 4 methods
%   4. measure hypergraph disimilarity
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 16, 2022

% function driver() %n, e, k)

% Set parameters
n = 20;
e = 40;
k = 3;

% 1. construct random hypergraph incidence matrix
W = randomIncidenceMatrix(n,e,k);
A = double(tensor(incidence2Adjacency(W)));

% H = Hypergraph();
% H.IM = W;

% 2. clique expand to generate graph data
D = cliqueExpand(W);

% 3. reconstruct hypergraph using 4 methods

% Compute multirelation tensor
multirelationMethods = ["zvi", "wang-zheng", "taylor", "jpic1"];
multirelationTensors = cell(4, 1);
for m=1:length(multirelationMethods)
    disp(multirelationMethods{m})
    T = multirelationTensor(D, k, multirelationMethods(m));
    multirelationTensors{m} = double(tensor(T));
end

% Map tensor to incidence matrix
%{
Wmatrices = cell(4, 1);
for m=1:length(multirelationMethods)
    w = zeros(n, nchoosek(n, k));
    idxs = nchoosek(1:n, k);
    T = multirelationTensors{m};
    for i=1:length(idxs)
        w(idxs(i,:), i) = T(idxs(i,:));
    end
    [~, ord] = sort(sum(w), 'descend');
    w = w(:,ord);
    w = (w > 1);
    Wmatrices{m} = w;
end
%}

% 4. measure hypergraph disimilarity
for m=1:length(multirelationMethods)

    % Hm = Hypergraph;
    % Hm.IM = Wmatrices{m};
    
    d = DissimilarityMeasures.hypergraphDissimilarity(A, multirelationTensors{m}, 'Hamming');
    disp(d);
    d = DissimilarityMeasures.hypergraphDissimilarity(A, multirelationTensors{m}, 'Spectral-H');
    disp(d);
    % d = DissimilarityMeasures.hypergraphDissimilarity(A, multirelationTensors{m}, 'Spectral-S');
    % disp(d);
    
end

% end
