function HDI = hubDepressed_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    HDI = sum(i_neighbors & j_neighbors) / max([sum(i_neighbors), sum(j_neighbors)]);
end
