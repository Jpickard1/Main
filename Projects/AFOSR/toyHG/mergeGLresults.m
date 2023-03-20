%% MERGE GREAT LAKES RESULTS
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 20, 2023

clear

tbls = cell(3,1);
load('toyHG/hyperchain.mat'); tbls{1} = r('hyperchain');
load('toyHG/hyperring.mat'); tbls{2} = r('hyperring');
load('toyHG/hyperstar.mat'); tbls{3} = r('hyperstar');

% Table of cells of strings to table of strings
% tbls2 = cell(3,1);
for k=1:3
    tt = tbls{k};
    ttt = table();
    for i=1:size(tt,1)
        for j=1:size(tt,2)
            val = tt{i,j};
            val = val{1};
            if isempty(val)
                val = "N/A";
            end
            ttt{i,j} = val
        end
    end
    ttt.Properties.VariableNames = tt.Properties.VariableNames
    tbls{k} = ttt;
end

expTable2latex(tbls{1}, 'null.tex')

