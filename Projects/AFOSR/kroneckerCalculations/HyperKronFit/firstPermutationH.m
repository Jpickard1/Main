function [p, P]=firstPermutationH(A, theta, maxItrs, verbose)
%FIRSTPERMUTATIONH
%
%   This sets the first permutation for the first iteration of the sample
%   gradient algorithm. Subsequent updates of the permutation are handled
%   through naxtPermutation.m
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

if nargin < 3
    maxItrs = 10000;
end
if nargin < 4
    verbose = false;
end

n = size(A,1);
p = randperm(n);

if nargout == 2
    P = p;
end

i = 0;
swap = 1;
while i < maxItrs
    j = randi([1 n]);
    k = randi([1 n]);

    u = rand();
    v2 = hPPRtest2(p, theta, A, j, k);
    if log(u) > v2
        if verbose
            disp(swap); swap = swap + 1;
            disp(p);
        end
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
        if nargout == 2
            P = [P; p];
        end
    end
    i = i + 1;
end
if verbose
    disp(p);
end

end
