%% Sample Permutation
%
%   This sets the first permutation for the first iteration of the sample
%   gradient algorithm. Subsequent updates of the permutation are handled
%   through naxtPermutation.m
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023
function [p]=firstPermutation(A, theta, maxItrs)

if nargin == 2
    maxItrs = 10000;
end

n = size(A,1);
p=1:n;

i = 0;
while i < maxItrs
    j = randi([1 n]);
    k = randi([1 n]);

    u = rand();
    v = PPRtest(p, theta, A, j, k);
    if u > v
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
    end
    i = i + 1;
end

end
