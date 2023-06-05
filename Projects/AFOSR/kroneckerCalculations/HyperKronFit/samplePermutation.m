%% Sample Permutation
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [p]=samplePermutation(theta,A)

n = size(A,1);

p = 1:n;
i = 1;
while conditionHere
    j = randi([1 n]);
    k = randi([1 n]);
    
    if rand() <  
        p([j k]) = p([k j]); % Swap elements j and k in the permutation
    end

end

end
