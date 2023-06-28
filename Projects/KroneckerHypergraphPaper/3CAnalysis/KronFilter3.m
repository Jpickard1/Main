function [IM, F] = KronFilter3(IM, levelSize)
%KRONFILTER3
%
% NORMALIZING IM, B, and C
%
%   Consider the following cases:
%       1. levelSize = size(IM,1):  Then IMp is 1x1 and size(Fp) = size(IM).
%       2. levelSize = 1: Then size(IMp) = size(IM) and Fp is 1x1.
%
%       In both cases, we must choose normalizations so that IM = IMp * Fp
%   
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 28, 2023

IMp = zeros(floor(size(IM,1) / levelSize), floor(size(IM,1) / levelSize));
for i=1:size(IMp)
    for j=1:size(IMp)
        IMp(i,j) = sum(IM((i-1)*levelSize+1:i*levelSize, (j-1)*levelSize+1:j*levelSize), 'all');
    end
end

% Compute F such that ||A - IMp kron F|| is minimized
Bvec = IMp(:);
n = levelSize;

% Create the kronecker product matrix for B
kronB = kron(eye(n*n), Bvec);

Ap = kron2vecPerm(IM);
Avec = Ap(:);

% Solve for C using lsqnonneg
C_vec = lsqnonneg(kronB, Avec);

% Reshape C to the appropriate size
Cnew = reshape(C_vec, [n n])';

F = {Cnew};
IM = {IMp};

end


