function HPI = hubPromoted_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    HPI = sum(i_neighbors & j_neighbors) / min([sum(i_neighbors), sum(j_neighbors)]);
end
