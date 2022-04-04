% Description: This function gives a popularity score for an edge based on
% equation 12 in: Uncovering missing links with cold ends by Zhu
function pop=edge_popularity(adj, i, j)
    if ~adj(i, j)
        pop = 0;
    else
        pop = (sum(adj(i,:)) - 1) * (sum(adj(j,:)) - 1);
    end
end
