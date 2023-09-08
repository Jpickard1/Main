function [A] = HG2Clique(HG)
%HG2Clique: constructs a clique graph from a hypergraph
%
%   HG = readEdgeSet('ES1ex.txt')
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

n = numel(HG.V);
A = sparse(n,n);
numE = size(HG.E,1);

for i=1:numE
    et = HG.E{i,1};
    eh = HG.E{i,2};
    if iscell(et)
        et = et{1};
    end
    if iscell(eh)
        eh = eh{1};
    end    
    for j=1:length(et)
        for k=1:length(eh)
            A(et(j),eh(k)) = 1;
        end
    end
end

end

