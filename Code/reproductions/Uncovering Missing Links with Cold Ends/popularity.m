% Description: This function gives a popularity score for an edge
function pop=popularity(adj, i, j)
    pop = (sum(adj(i,:)) - 1) * (sum(adj(j,:)) - 1);
end