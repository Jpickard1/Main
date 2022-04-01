function RA = resourceAllocation_index(adj, i, j)
    i_neighbors = (adj(i,:) == 1);
    j_neighbors = (adj(j,:) == 1);
    union = i_neighbors & j_neighbors;
    RA = 0;
    for k=1:length(adj)
        if union(k)
            RA = RA + (1/sum(adj(k,:)));
        end
    end
end
