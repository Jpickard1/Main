function [A, removed_vxs] = reduceConnComp(A)
%REDUCECONNCOMP This function reduces a graph to the single largest 
%   connected component
%
% Auth: Joshua Pickard (jpic@umich.edu)
% Date: May 23, 2022
bins = conncomp(graph(A));
largest_component = mode(bins);

keep = bins == (largest_component);
A = A(keep, keep);

removed_vxs = find(bins ~= largest_component);
end

