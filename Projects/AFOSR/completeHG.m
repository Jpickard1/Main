function [HG] = completeHG(V, k)
%HYPERRING This function returns a complete hypergraph
%
%   A copmlete (hyper)graph contains all possible edges.
%
% PARAMETERS:
%   V: number of vertices
%   k: order of hypergraph
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 11, 2023

E = nchoosek(1:V,k);
IM = HAT.hyperedges2IM(E);
HG = Hypergraph('IM', sparse(IM));

end

