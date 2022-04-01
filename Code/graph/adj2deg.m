% Description: This function takes an adjacency matrix and produces the
% corresponding degree matrix of the graph.
function deg = adj2deg(adj)
    deg = diag(sum(adj));
end