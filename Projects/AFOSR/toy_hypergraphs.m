%% Toy Hypergraph
%
%   This file computes the minimal set of observable nodes for a set of
%   small hypergraphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 2, 2023
%
% MOD: March 16, 2023 - The initial draft performed the experiment, but was
%                       not thoughtfullying saving results to a file or
%                       table. Additionally, only the number of MON and not
%                       actual MON were saved. Updates were made to save
%                       the specific set of MON to a file.
%                     - Functionized this code to run as array job on great
%                       lakes

function toy_hypergraphs(tiArr)

% Set of possible parameters
N=3:6;
K=3:5;
type = ["hyperring", "hyperchain", "hyperstar"];

r = containers.Map;
for ti=1:length(type)
    t = type(ti);
    r(t) = table(); % rows are number of vertices, columns are order k
end

tic;
for ki=1:length(K)
    k = K(ki);
    for ni=1:length(N)
        n = N(ni);
        if n < k
            continue
        end
        % for ti=1:length(type)
            t = type(tiArr);

            % Print parameters
            disp(string(t) + "(" + string(n) + "," + string(k) + ")");

            % Experiment
            HG = getToyHG(n,k,t);               % Get hypergraph
            O = getObservabilityMatrices(HG);   % Compute individual observability matrices
            [D, ~] = greedyMON(O, n);             % Greedy Node Selection

            T = r(t);
            T{ni, ki} = {D};
            r(t) = T;

            % Save result
            disp(D);

            fileName = "toyHG/" + string(t) + ".mat";
            cmd = "save " + fileName + " " + "r -v7.3";
            eval(cmd);
            disp(cmd);

        % end
    end
end
disp(toc);

end

