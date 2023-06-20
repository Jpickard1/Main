function [E, F] = vrcomplex(points, r)
% Compute the Vietoris-Rips complex for a given point cloud using 
% Boissonnat and Maria's 2013 algorithm [1].
% 
% Inputs:
% - points: a n x d matrix representing n points in d-dimensional space
% - r: the radius parameter for the Vietoris-Rips complex
%
% Outputs:
% - E: a m x 2 matrix representing the edges (m is the number of edges)
% - F: a p x 3 matrix representing the triangles (p is the number of triangles)
%
% References:
% - [1] Computing Persistent Homology with Various Coefficient Fields in a 
%       Single Pass
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 12, 2023

n = size(points, 1);
d = size(points, 2);

% Create a kd-tree for the point cloud
tree = KDTreeSearcher(points);

% Initialize the output variables
E = zeros(0, 2);
F = zeros(0, 3);

% Process each point in the point cloud
for i = 1:n
    % Find all points within radius r of point i
    idx = rangesearch(tree, points(i,:), r);
    idx = idx{1};
    
    % Add edges between point i and its neighbors
    new_edges = [repmat(i, length(idx), 1), idx(:)];
    E = [E; new_edges];
    
    % Check if new triangles can be formed
    for j = 1:length(idx)-1
        for k = j+1:length(idx)
            if norm(points(idx(j),:)-points(idx(k),:)) <= 2*r
                % Add triangle between point i and its neighbors
                new_face = [i, idx(j), idx(k)];
                F = [F; new_face];
            end
        end
    end
end

% Remove duplicate edges and triangles
E = unique(sort(E, 2), 'rows');
F = unique(sort(F, 2), 'rows');
end
