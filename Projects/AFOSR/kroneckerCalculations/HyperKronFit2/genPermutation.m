function [p, P, idxs]=genPermutation(p, theta, maxItrs, E)
%GENPERMUTATION
%
%   This sets the first permutation for the first iteration of the sample
%   gradient algorithm. Subsequent updates of the permutation are handled
%   through nextPermutation.m
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

n = max(max(E));
% n = size(A,1);

if nargout == 2
    P = zeros(round(maxItrs * 0.05), n);
    P(1,:) = p;
end

i = 0;
swap = 0;
while i < maxItrs
    j = randi([1 n]);
    k = randi([1 n]);
    while j == k; k = randi([1 n]); end
    idx = [j k];

    v = permutationProbabilityRatioTest(p, idx, E, theta);

    u = rand();
    if log(u) > v
        swap = swap + 1;
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
        if nargout == 2; P(swap,:) = p; end
    end
    i = i + 1;
end

if nargout == 2
    P = P(1:swap,:);
end
if nargout == 3
    P = (swap > 0);
    idxs = [j k];
end
if maxItrs == 1
    P = (swap > 0);
end

end
