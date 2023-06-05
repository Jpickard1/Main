function [r] = permutationProbabilityRatio(p1, p2, theta, A)
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

n = size(A,1);
k = log(n) / log(2);

r = 1;
for i=1:n
    for j=1:n
        e1 = edgeProbability(n, theta, p1(i), p1(j));
        e2 = edgeProbability(n, theta, p2(i), p2(j));
        if A(i,j) == 1
            rij = e1 / e2;
        else
            rij = (1 - e1) / (1 - e2);
        end
        r = r * rij;
    end
end

end