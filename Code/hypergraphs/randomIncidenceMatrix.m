%% Multirelations
%   This function constructs a random incidence matrix with n vertices, e
%   edges, and k vertices per edge.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 16, 2022
function W = randomIncidenceMatrix(n,e,k)
    W = zeros(n, e);
    for i=1:e
        vxc = randperm(n,k);
        W(vxc, i) = 1;
    end
end
