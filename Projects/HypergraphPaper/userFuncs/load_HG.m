function HG = load_HG(ds)
%LOAD_HG Summary of this function goes here
%   Detailed explanation goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 13, 2023

% DS = {"DBLP", "Cora Reference", "Cora Citation", "Citeseer Reference", "Citeseer Citation"}
DS = {'DBLP', 'Cora Reference', 'Cora Citation', 'Citeseer Reference', 'Citeseer Citation', 'ArnetMiner Citation', 'Oddysey'};

if ~any(strcmp(DS,ds))
    error('INVALID DATA SET: ' + ds)
end

S = HAT.load(ds);
IM = full(S);
HG = Hypergraph('IM', IM);

end

