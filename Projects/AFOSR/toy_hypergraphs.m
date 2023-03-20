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
% %% Debug
% tiArr = 3;
maxItr = 100;

addpath(genpath('/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/'))
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/'))

% Set of possible parameters
N=3:8; colNames = cell(length(N), 1); for i=1:length(N); colNames{i} = char(string(N(i))); end
K=2:7;
type = ["hyperring", "hyperchain", "hyperstar"];

r = containers.Map;
for ti=1:length(type)
    t = type(ti);
    r(t) = cell2table(cell(length(K), length(N)), 'VariableNames', colNames); % rows are number of vertices, columns are order k
    % r(t).Properties.RowNames(1:length(N))
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
            HG = getToyHG(n,k,t)               % Get hypergraph
            for itr=1:maxItr
                x = rand(n,1);
                O = HGObsvNum(HG, x);
                [d, ~] = greedyMON(O, n);
                dStr = strjoin(string(d), ' ');
                if itr==1; D = dStr;
                elseif ~contains(D, dStr)
                    D = D + "/" + dStr;
                end
            end
            % O = HGObsvNum(HG, rand(n,1));
            % O = getObservabilityMatricesNumeric(HG, rand(n,1));   % Compute individual observability matrices
            % [D, ~] = greedyMON(O, n);             % Greedy Node Selection

            T = r(t);
            T{ki, ni} = {D};
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

