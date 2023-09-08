%% Code Assist
%   This code was generated with Chat GPT and the Python code published
%   with the paper.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

objective_manyIC(rand(3,1), rand(3,3), rand(3,1))

%% Functions

function obj = objective_manyIC(x, At, x0)
    % The objective is to maximize the output energy
    % x is the design variable or the output gain
    % At is the KO raised to the power t
    % x0 is the initial condition
    
    obj = 0;
    for i = 1:numel(At)
        obj = obj - norm(x * At{i} * x0, 2); % for use with energy_max or energy_max_uni_dist
        % obj = obj - norm(reshape(x * At{i} * x0, 1, 1), 2); % for use with energy_max_with_mean
    end
end

function [C, r] = energy_maximization_with_uni_dist(X, A, ntimepts, Tf, IC)
    disp('------Optimizing for gene sampling weights------');
    start_time = tic;
    r = randi(100);
    rng(r);
    C0 = normrnd(4000, 1000.0, 1, length(A)); % initialization
    At = cell(1, Tf + 1);
    for i = 0:Tf
        At{i + 1} = A^i;
    end
    x0 = X(:, 1, :);
    % get the min and max of each gene's IC
    x0min = min(x0, [], 2);
    x0max = max(x0, [], 2);
    % form a set of IC's distributed uniformly
    numICs = size(X, 3);
    x0uni = zeros(length(x0min), numICs);
    for ii = 1:numICs
        x0tmp = unifrnd(x0min, x0max);
        x0uni(:, ii) = x0tmp;
    end
    disp(['Initial objective: ', num2str(objective_manyIC(C0, At, x0uni))]);
    % optimize
    lb = zeros(1, size(X, 1)); % C should be nonnegative
    ub = [];
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
    solution = fmincon(@(C) objective_manyIC(C, At, x0uni), C0, [], [], [], [], lb, ub, [], options);
    C = solution;
    % show final objective
    disp(['Final objective: ', num2str(objective_manyIC(C, At, x0uni))]);
    disp([num2str(toc(start_time) / 60), ' minutes']);
end

function rho = reconstruct_x0_sequential(data_fc_norm, nT, A, C, CsortedInds, samplingFreq, order)
    % nT is the number of timepoints for which to generate outputs and reconstruct x0 with.
    % First sample genes by rank (or randomly), compute reconstruct error, 
    % then sample again (adding to the already sampled genes)
    % and repeat
    
    topNums = 1:samplingFreq:length(data_fc_norm); % sampling these genes in order of rank samplingFreq at a time
    % don't need to sample all genes, as we know the reconstruction loss asymptotes after 80 genes or so.
    rho = [];
    Uset = 1:length(data_fc_norm); % set of all genes
    NUset = Uset; % updated each iteration to remove the genes that were used prior (no double dipping)
    
    for ii = 1:(length(topNums) - 1)
        if strcmp(order, 'top')
            top_inds = CsortedInds(end:-1:(topNums(1):topNums(ii+1)));
            top_inds_complement = setdiff(1:length(data_fc_norm), top_inds);
        elseif strcmp(order, 'random')
            if length(NUset) <= samplingFreq
                break;
            end
            if ii == 1
                rs = datasample(NUset, samplingFreq, 'Replace', false);
                store_rs = rs;
            else
                rs = datasample(NUset, samplingFreq, 'Replace', false);
                store_rs = [store_rs, rs];
            end
            NUset = setdiff(NUset, store_rs);
            top_inds_complement = NUset;
        end
        
        C_sensors = C;
        C_sensors(:, top_inds_complement) = 0.0;
        
        % generate outputs at nT timepoints, check rank of O_T and estimate IC
        O_T = zeros(length(C_sensors) * nT, size(C_sensors, 2)); % O_T has dim nOutputs*nTimepoints,nStateVars
        for jj = 1:size(O_T, 1) % observability matrix
            O_T(jj, :) = C_sensors * (A^(jj - 1))';
        end

        % use the learned dynamics and sampling to generate output over time
        y_r1 = zeros(length(C), nT);
        y_r2 = zeros(length(C), nT);
        for jj = 1:nT
            y_r1(:, jj) = C_sensors * (A^(jj - 1))' * data_fc_norm(:, 1, 1);
            y_r2(:, jj) = C_sensors * (A^(jj - 1))' * data_fc_norm(:, 1, 2);
        end

        % estimate of x0 
        x0_est_r1 = pinv(O_T) * y_r1';
        x0_est_r2 = pinv(O_T) * y_r2';

        rho1 = corr(data_fc_norm(:, 1, 1), x0_est_r1', 'Type', 'Pearson');
        rho2 = corr(data_fc_norm(:, 1, 2), x0_est_r2', 'Type', 'Pearson');
        rho(end+1) = (rho1 + rho2) / 2;
    end
end

function rho = reconstruct_x0(data, nT, A, C)
    % The matrix C should be of shape (pxn) if the data are shape (nxm) for p measured genes, n total genes, and m samples
    % If p << n, leave only entries of C that correspond to the p genes as nonzero.

    % generate outputs at nT timepoints, check rank of O_T, and estimate IC
    O_T = zeros(length(C) * nT, size(C, 2)); % O_T has dim p*nT x n
    for ii = 1:size(O_T, 1) % observability matrix
        O_T(ii, :) = C * (A^(ii - 1))';
    end

    % use the learned dynamics and sampling to generate output over time
    y_r1 = zeros(length(C), nT);
    y_r2 = zeros(length(C), nT);
    for ii = 1:nT
        y_r1(:, ii) = C * (A^(ii - 1))' * data(:, 1, 1);
        y_r2(:, ii) = C * (A^(ii - 1))' * data(:, 1, 2);
    end

    % estimate of x0 
    x0_est_r1 = pinv(O_T) * y_r1';
    x0_est_r2 = pinv(O_T) * y_r2';

    rho1 = corr(data(:, 1, 1), x0_est_r1', 'Type', 'Pearson');
    rho2 = corr(data(:, 1, 2), x0_est_r2', 'Type', 'Pearson');
    rho = (rho1 + rho2) / 2;
end

