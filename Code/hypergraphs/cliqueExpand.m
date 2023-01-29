%% CLIQUE EXPAND
%   This function clique expands the incidence matrix of a hypergraph to 
%   construct the adjacency matrix of a graph.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 16, 2022
function A = cliqueExpand(W)

% [n, e] = size(W);
% A = zeros(n, n);
A = W * W';
A = (A > 0);
A = A - diag(diag(A));

end
