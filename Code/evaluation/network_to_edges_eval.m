% Description: This function accepts an adjacency matrix (network) and
% vectorizes it into a binary array (1 where there is an edge and 0
% otherwise). This is deterministic, to keep it consistent, and this
% process only considers each edge once, by only looking at the edges in
% the upper triangular part of the matrix (above the main diagonal). This
% function is intended to be used in eval_metrics, to vectorize the
% predicted links and compute performance metrics.
%
% Created: Joshua Pickard, 4/1/2022
%
function edges = network_to_edges_eval(network)
    edges = [];
    [m, n] = size(network);
    if m ~= n
        disp('WARNING: Network matrix is not symmetric')
    end
    for i=1:m
        for j=i+1:n
            edges = [edges; network(i,j)];
        end
    end
end
