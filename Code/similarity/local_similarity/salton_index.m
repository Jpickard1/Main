% Description: Returns the salton similarity index between nodes i and j 
% based on the adjacency matrix. The salton index is the cardinality of the
% union between the set of neighbors of the two nodes divided by the square
% root of the product of the degree of each node.
%
% Salton: |neighbors_i U neighbors_j| / sqrt(degree_i * degree_j)
%
% Parameters:  adj is the adjacency matrix. It must be a logical matrix,
% where 1 indicates a link and 0 indicates no link.
function salton=salton_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    salton = sum(i_neighbors & j_neighbors) / sqrt(sum(i_neighbors) * sum(j_neighbors));
end
