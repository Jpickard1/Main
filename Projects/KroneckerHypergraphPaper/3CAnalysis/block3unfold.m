function [M] = block3unfold(A, b)
%BLOCK3UNFOLD Unfolds 3-way tensor using tensor blocks
%
%   I don't know if this follows exactly the method proposed by Regarisson
%   and Van Loan, but you should probably read their paper Tensor Block
%   Unfoldings
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 29, 2023

n = size(A, 1);
M = zeros(n, n^2);

for i=1:(n/b)
    for j=1:(n/b)
        for k=1:(n/b)
            unfoldedBlock = reshape(A((i-1)*b+1:i*b, (j-1)*b+1:j*b, (k-1)*b+1:k*b), [b,b^2]);
            M((i-1)*b+1:i*b, (j-1)*b^3+(k-1)*b^2+1:(j-1)*b^3+k*b^2) = unfoldedBlock;
        end
    end
end

end

