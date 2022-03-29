function AA = adamicAdar_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    union = i_neighbors & j_neighbors;
    AA = 0;
    for k=1:length(adj)
        if union(k)
            AA = AA + (1/log(sum(adj(k,:))));
        end
    end
end
