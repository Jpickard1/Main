function [A] = E2A(V, E)
%E2A Constructe adjacency matrix from edge set
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 9, 2023

A = zeros(V,V);
for i=1:size(E,1)
    A(E(i,1), E(i,2)) = 1;
end
A = A + A';
A = logical(A);

end

