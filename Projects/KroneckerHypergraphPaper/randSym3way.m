function [A] = randSym3way(n)
%RANDSYM3WAY Returns a random 3 way symmetric tensor
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 1, 2023

A = zeros(n,n,n);
for i=1:n
    for j=i:n
        for k=j:n
            r = rand();
            A(i,j,k) = r;
            A(i,k,j) = r;
            A(j,i,k) = r;
            A(j,k,i) = r;
            A(k,i,j) = r;
            A(k,j,i) = r;
        end
    end
end

end