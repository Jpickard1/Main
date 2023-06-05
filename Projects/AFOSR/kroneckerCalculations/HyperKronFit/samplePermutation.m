%% Sample Permutation
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [p]=samplePermutation(A,theta)

n = size(A,1);

p = 1:n;
c = 0;
i = 0;
while c < 100 && i < 10000
    j = randi([1 n]);
    k = randi([1 n]);

    u = rand();
    v = permutationProbabilityRatio(p, theta, A, j, k);

    if u < v
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
        i = i + c;
        c = 0;
    end
    c = c + 1;
end

end
