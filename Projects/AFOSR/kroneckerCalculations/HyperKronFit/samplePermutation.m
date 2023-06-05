%% Sample Permutation
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [p]=samplePermutation(A,theta)

n = size(A,1);

p = 1:n;
i = 1;
while i < 10
    j = randi([1 n]);
    k = randi([1 n]);

    u = rand();
    p2 = p; p2([j k]) = p2([k j]);
    v = permutationProbabilityRatio(p, p2, theta, A);

    if u < v
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
    end
    i = i + 1;
end

end
