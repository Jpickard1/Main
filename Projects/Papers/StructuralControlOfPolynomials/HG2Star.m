function [A] = HG2Star(HG)
%HG2Star: constructs a star graph from a hypergraph
%
%   HG = readEdgeSet('ES1ex.txt')
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

n = numel(HG.V);
numE = size(HG.E,1);
A = sparse(n+numE,n+numE);

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
        A(et(j),n+i) = 1;
    end
    for j=1:length(eh)
        A(n+i,eh(j)) = 1;
    end
end
full(A)

end