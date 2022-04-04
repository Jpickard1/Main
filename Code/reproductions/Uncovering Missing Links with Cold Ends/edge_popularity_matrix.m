function matrix=edge_popularity_matrix(adj)
    matrix = zeros(size(adj));
    for i=1:length(adj)
        for j=i+1:length(adj)
            pop = edge_popularity(adj, i, j);
            matrix(i, j) = pop;
            matrix(j, i) = pop;
        end
    end
end