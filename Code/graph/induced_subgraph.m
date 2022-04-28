function subgraph = induced_subgraph(A, vxs)
% SUBGRAPH: This function returns an induced subgraph on the vertices
% specified in vxs
    subgraph = A(:, vxs);
    subgraph = subgraph(vxs, :);
end

%{
Sample parameters to use
A = randi([0, 1], 5,5);
A = A + A';
A = logical(A);
vxs = [1 3 5];
%}