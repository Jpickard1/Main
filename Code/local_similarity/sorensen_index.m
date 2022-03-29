function sorensen=sorensen_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    sorensen = sum(i_neighbors & j_neighbors) / (sum(i_neighbors) + sum(j_neighbors));
end
