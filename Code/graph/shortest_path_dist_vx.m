function [distances] = shortest_path_dist_vx(A, vx, known)
%SHORTEST_PATH_DIST_VX This function computes the minimal distance from a
%   vx to all other vxs in a graph
%
% Auth: Joshua Pickard
% Date: May 10, 2022

% Save distance data
distances = zeros(length(A),1);

% Record which nodes have been found
found = false(length(A),1);
found(vx) = true;
if nargin == 3
    found = found | known;
end

% Record frontier nodes
frontier = zeros(length(A), 1);
frontier(vx) = 1;

len = 1;
prev_sum = 0;
while sum(found) ~= prev_sum
    prev_sum = sum(found);
    frontier_new = A * frontier;
    new_nodes = logical(frontier_new) & ~found;
    found(new_nodes) = true;
    distances(new_nodes) = len;
    len = len + 1;
    frontier = new_nodes;
end

% Reset this value
distances(vx) = 0;
end

