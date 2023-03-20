%% MERGE GREAT LAKES RESULTS
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 20, 2023

%% merge numeric results
clear
tbls = cell(3,1);
load('hyperchain.mat'); tbls{1} = r('hyperchain');
load('hyperring.mat'); tbls{2} = r('hyperring');
load('hyperstar.mat'); tbls{3} = r('hyperstar');

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

expTable2latex(tbls{1}, 'hyperchain.tex')
expTable2latex(tbls{2}, 'hyperring.tex')
expTable2latex(tbls{3}, 'hyperstar.tex')

%% symbolic results
clear
tbls = cell(3,1);
load('hyperchain_sym.mat'); tbls{1} = r('hyperchain');
load('hyperring_sym.mat'); tbls{2} = r('hyperring');
load('hyperstar_sym.mat'); tbls{3} = r('hyperstar');

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