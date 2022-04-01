function PA = preferentialAttachment_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    PA = (sum(i_neighbors) * sum(j_neighbors));
end
