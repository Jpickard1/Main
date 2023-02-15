function [scores, labels] = predictHyperedges(HGtrue, HGobserved, F)
%PREDICTHYPEREDGES 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: January 29, 2023

%{
HGtrue = HGt;
HGobserved = O1{6};
%}

%% TODO: set this section to run as a function
%% Params
nPredictors = 7;    % Number of predicotrs to use
predictors = ["CN", "LHN", "SA", "SO", "HP", "HD", "JA"];   % List of predictors

% Extract incidence matrices of true and observed hypergraphs
IMt = full(HGtrue.IM);
IMo = full(HGobserved.IM);

% Get indices of known and unknown hyperedges in the incidence matrices
Kidx = find(sum(IMo, 1) >= 0.5 * sum(IMt, 1));          % Get idxs of Known hyperedges
Uidx = setdiff(1:size(IMt, 2), Kidx);  % Get idxs of Unknown hyperedges

% Get unknown (True) and False hyperedges for prediction
U = trueHyperedges(HGtrue, Uidx);   % True hyperedge cell array
% F = falseHyperedges(HGtrue, Uidx);  % False hyperedge cell array
H = [F; U];                             % True + False hyperedges to evaluate

labels = ones(length(H), 1);            % Set True labels
labels(1:length(F)) = 0;                % Set False labels
scores = zeros(length(H), nPredictors); % Scores of each hyperedge for each predictor

%% Compute pairwise similarity matrices betweeen vertices

% Construct unweighted clique graph
C = IMo * IMo';
% C(C > 0) = 1;

similarity = cell(nPredictors, 1);
for i=1:nPredictors
    % disp(predictors(i));
    if strcmp(predictors(i), "CN")
        similarity{i} = C * C;
    elseif strcmp(predictors(i), "LHN")
        CN = similarity{1};
        d = sum(C);
        D = d' * d;
        similarity{i} = CN ./ D;
    elseif strcmp(predictors(i), "SA")
        CN = similarity{1};
        d = sum(C);
        D = d' * d;
        D = sqrt(D);
        similarity{i} = CN ./ D;
    elseif strcmp(predictors(i), "SO")
        CN = similarity{1};
        d = sum(C);
        D = d' + d;
        similarity{i} = 2 * CN ./ D;
    elseif strcmp(predictors(i), "HD")
        CN = similarity{1};
        d = sum(C);
        r1 = repmat(d,length(d), 1);
        r2 = repmat(d',1,length(d));
        D = max(r1, r2);
        similarity{i} = 2 * CN ./ D;
    elseif strcmp(predictors(i), "HP")
        CN = similarity{1};
        d = sum(C);
        r1 = repmat(d,length(d), 1);
        r2 = repmat(d',1,length(d));
        D = min(r1, r2);
        similarity{i} = 2 * CN ./ D;
    elseif strcmp(predictors(i), "JA")
        CN = similarity{1};
        TN = zeros(size(CN));
        N = cell(length(CN), 1);
        for j=1:length(CN)
            N{j} = find(C(i,:) ~= 0);
        end
        for j=1:length(CN)
            for k=j+1:length(CN)
                TN(j,k) = length(find(N{j} + N{k} ~= 0));
%                TN(j,k) = length(find(((C(i,:) ~= 0) + (C(j,:) ~= 0)) ~= 0));
                TN(k,j) = TN(j,k);
            end
        end
        similarity{i} = CN ./ TN;
    end
    similarity{i}(isnan(similarity{i})) = 0;
    similarity{i}(isinf(similarity{i})) = 0;
    % similarity{i} = compute_similarity(C, predictors(i));
end

%% Score Hyperedges

ss = size(similarity{1});

% parfor i=1:length(H)
for i=1:length(H)
    e = H{i};
    p = nchoosek(e, 2);
    vecP = sub2ind(ss, p(:,1), p(:,2));
    % disp(length(vecP));
    % Loop over all edge predictors
    for j=1:nPredictors
        scores(i,j) = sum(similarity{j}(vecP)) / length(vecP);        
        % disp(string(i) + ":" + string(j));
    end
end


end

