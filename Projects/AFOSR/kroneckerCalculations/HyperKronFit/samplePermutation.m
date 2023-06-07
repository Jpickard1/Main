%% Sample Permutation
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [p]=samplePermutation(A,theta,p)

n = size(A,1);

if nargin == 2
    p = 1:n;
end

c = 0;
i = 0;
while c < 10 && i < 100
    j = randi([1 n]);
    k = randi([1 n]);
    while j==k; k = randi([1 n]); end 

    u = rand();
    v = PPRtest(p, theta, A, j, k);
    % v = permutationProbabilityRatio(p, theta, A, j, k);

    if u > v
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
        i = i + c;
        c = 0;
    end
    c = c + 1;
end

end
