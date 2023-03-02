%% Toy Hypergraph
%
%   This file computes the minimal set of observable nodes for a set of
%   small hypergraphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 2, 2023

% Set of possible parameters
N=3:7;
K=3:5;
type = ["hyperring", "hyperchain", "hyperstar"];

tic;
for ki=1:length(K)
    k = K(ki);
    for ni=1:length(N)
        n = N(ni);
        for ti=1:length(type)
            t = type(ti);

            % Print parameters
            disp(string(t) + "(" + string(n) + "," + string(k) + ")");

            % Experiment
            HG = getToyHG(n,k,t);               % Get hypergraph
            O = getObservabilityMatrices(HG);   % Compute individual observability matrices
            [D, OD] = greedyMON(O, n);             % Greedy Node Selection

            % Save result
            disp(D);

        end
    end
end
disp(toc);

