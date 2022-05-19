function gamma = clustering_coef_vx(A, vx, loops)
% CLUSTERING COEF VX: This function returns the clustering coefficient for
%   single vetrex (vx) in a graph. In Small Worlds (Watts 1998), the
%   clustering coefficient is defined
%
%       gamme = (# of edges between neighbors) / (# possible edges)
%
%   NOTES
%       - loops: if the graph is a simple graph, the loops should be set to
%       false. This does not allow for self loops. Otherwise set loops to
%       true.
%
% Auth: Joshua Pickard
% Date: May 10, 2022

    % Select neighbors
    neighbors = find(A(vx,:));
    if sum(ismember(neighbors, vx)) == 0
        neighbors = [neighbors vx];
    end
    % Extract induced subgraph
    net = A(neighbors,:);
    net = net(:, neighbors);

    % Compute number of edges
    num_edges = ((sum(sum(net)) - sum(diag(net))) / 2) + sum(diag(net));
    if length(neighbors) < 2
        disp("WARNING: in clustering_coef_vx: length(neighbors) < 2");
        disp("         gamma = 0 returned");
        disp("         vx = " + string(vx));
        gamma = 0;
        return;
    end
    possible_edges = nchoosek(length(neighbors), 2);
    if loops
        possible_edges = possible_edges + length(neighbors);
    end
    % Compute clustering coefficient
    gamma = num_edges / possible_edges;
end
