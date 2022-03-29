% Description: Returns the common neighbors similarity index (CN) between
% nodes i and j based on the adjacency matrix. The CN index is the 
% cardinality of the union between the set of neighbors of the two nodes.
%
% Parameters:  adj is the adjacency matrix. It must be a logical matrix,
% where 1 indicates a link and 0 indicates no link.
function CN=commonNeighbors(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    CN = sum(i_neighbors & j_neighbors);
end
