function CPL = characteristic_path_length(A)
% CHARACTERISTIC_PATH_LENGTH This function computes the characteristic path
%   length of the graph defined with the adjacency matrix A.
%
% Auth: Joshua Pickard
% Date: May 12, 2022

%%
shortest_path = zeros(length(A));
known = false(length(A), 1);

for vx=1:length(A)
    disp(string(vx));
    shortest_path(vx,:) = shortest_path_dist_vx(A, vx, known);
    known(vx) = true;
end

shortest_path = shortest_path + shortest_path';

% Compute characteristic path length
avg_vx_dist = mean(shortest_path);
CPL = median(avg_vx_dist);

end
