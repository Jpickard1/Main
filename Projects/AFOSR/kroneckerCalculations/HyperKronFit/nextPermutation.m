function [p,accept, j, k]=nextPermutation(A,theta,p)
%NEXTPERMUTATION
%                   
%   This generates the next permutation used by sampleGradient.m from a
%   prior permutation
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

n = size(A, 1);

% Select 2 random indices to change
j = randi([1 n]);
k = randi([1 n]);

% Generate random number to regulate the acceptance of the next permutation
u = rand();

% Calculate the likelihood ratio
v = PPRtest2(p, theta, A, j, k);

% Check if the updated permutation is accepted
accept = false;
if log(u) > v
    p([j k]) = p([k j]); % Swap elements j and k in the permutation
    accept = true;
end

end
