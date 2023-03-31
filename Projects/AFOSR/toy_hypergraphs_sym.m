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
% MOD: March 20, 2023 - I updated the methed of computing the derivatives
%                       I made this file evaluate the derivatives
%                       numerically, up to 100 iterations. I made a new
%                       file, toy_hypergraphs_sym.m that performs the
%                       symbolic calculations.
% MOD: March 22, 2023 - The symbolic computations for n=8 are going too
%                       slow to have nice results by Friday, so I reduced
%                       the bounds to n=6, changed the output file name,
%                       and I am rerunning this job

function toy_hypergraphs_sym(tiArr)
% %% Debug
% tiArr = 3;
maxItr = 100;

addpath(genpath('/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/'))
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/'))

% Set of possible parameters
N=3:6; colNames = cell(length(N), 1); for i=1:length(N); colNames{i} = char(string(N(i))); end
K=5:7;
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
            O = HGObsvSym(HG);
            % O = getObservabilityMatricesNumeric(HG, rand(n,1));   % Compute individual observability matrices
            [D, ~] = greedyMON(O, n);             % Greedy Node Selection

            T = r(t);
            T{ki, ni} = {D};
            r(t) = T;

            % Save result
            disp(D);

            fileName = "toyHG/" + string(t) + "_sym_2.mat";
            cmd = "save " + fileName + " " + "r -v7.3";
            eval(cmd);
            disp(cmd);

        % end
    end
end
disp(toc);

end

