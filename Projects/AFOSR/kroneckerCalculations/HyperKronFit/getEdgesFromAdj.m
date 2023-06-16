function [E] = getEdgesFromAdj(A)
%GETEDGESFROMADJ
%
%   This function gets a list of all edges in a graph and hyperedges in a
%   hypergraph from an adjacency tensor (or matrix). If k=2 and it is an
%   adjacency matrix, then the edge list is directed. If k>2 and it is an
%   adjacency tensor for a hypergraph, then the edges are undirected and no
%   duplicate edges are returned.
%
%   Particularly for hypergraphs, A must be sparse with the ndSparse class
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 15, 2023

if ndims(A) == 2
    [e1, e2] = find(A ~= 0);
    E = [e1 e2];
else
    E = getHyperedges(A);
end

end
