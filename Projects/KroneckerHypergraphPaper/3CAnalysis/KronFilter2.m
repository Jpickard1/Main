function [IM, F] = KronFilter2(IM, levelSize)
%KRONFILTER
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
    Fp  = zeros(levelSize, levelSize);
    % Fp  = zeros(levelSize, levelSize, levelSize);
    for i=1:size(IMp)
        for j=1:size(IMp)
            IMp(i,j) = sum(IM((i-1)*levelSize+1:i*levelSize, (j-1)*levelSize+1:j*levelSize), 'all');
            if IMp(i,j) ~= 0
                Fp = Fp + (IM((i-1)*levelSize+1:i*levelSize, (j-1)*levelSize+1:j*levelSize) / (IMp(i,j)));
            else
                disp('skip')
            end
            % Fp(:,:,i) = IM((i-1)*levelSize+1:i*levelSize, (j-1)*levelSize+1:j*levelSize);
        end
    end
    % F = median(Fp,3);
    F = Fp ./ (i*j);
    IM = IMp;
end

