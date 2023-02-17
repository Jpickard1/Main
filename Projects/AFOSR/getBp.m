function B = getBp(Amat, p, k)
%GETBP
%
%   INPUTS:
%       Amat: unfolded adjacency tensor
%       p:    integer
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023

n = size(Amat, 1);

B = sparse(n^(p+1), size(Amat, 2) * n^(p+1));    % Set sparse matrix to return
bound = (p-1)*k-(2*p-3);        % Compute upper bounds of loop
for i=1:bound
    if i==1
        P = sparse(Amat);
    else
        P = sparse(eye(n,n));
    end
    for j=2:bound
        if j ~= i
            P = kron(P, eye(n,n));
        else
            P = kron(P, Amat);
        end
    end
    if i~=1
        B = B + P;
    else
        B = P;
    end
end

end
