function [r] = hPPRtest(p, theta, A, u, v)
% HPPRTEST
%
%   This function evaluates a swap of nodes u and v in the node permutation
%   p based on the Kronecker initiator theta and hypergraph adjacency
%   tensor A.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

n = size(A,1);          % number of vertices in A
K = length(size(A));    % order of uniform hypergraph
p1 = p;                 % permutation 1
p2 = p;
p2([u v]) = p2([v u]);  % permutation 2

% Iterate over modes
for mode=1:K

    % Iterate over vertices
    for vx=1:n

    end

end

end

% r = 1;
% for i=1:n
%     e1 = edgeProbability(n, theta, p1(u), p1(i));
%     e2 = edgeProbability(n, theta, p2(u), p2(i));
%     if A(u,i) == 1
%         rij = e1 / e2;
%     else
%         rij = (1 - e1) / (1 - e2);
%     end
%     r = r * (rij);
% 
%     e1 = edgeProbability(n, theta, p1(i), p1(u));
%     e2 = edgeProbability(n, theta, p2(i), p2(u));
%     if A(i,u) == 1
%         rij = e1 / e2;
%     else
%         rij = (1 - e1) / (1 - e2);
%     end
%     r = r * (rij);
% 
%     e1 = edgeProbability(n, theta, p1(v), p1(i));
%     e2 = edgeProbability(n, theta, p2(v), p2(i));
%     if A(v,i) == 1
%         rij = e1 / e2;
%     else
%         rij = (1 - e1) / (1 - e2);
%     end
%     r = r * (rij);
% 
%     e1 = edgeProbability(n, theta, p1(i), p1(v));
%     e2 = edgeProbability(n, theta, p2(i), p2(v));
%     if A(i,v) == 1
%         rij = e1 / e2;
%     else
%         rij = (1 - e1) / (1 - e2);
%     end
%     r = r * (rij);
% end