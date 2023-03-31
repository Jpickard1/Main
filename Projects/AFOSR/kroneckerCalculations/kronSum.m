function [s] = kronSum(x, y, n)
%KRONSUM 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 31, 2023

for i=1:n
    if i==1
        r = y;
    else
        r = x;
    end
    for j=2:n
        if i ~= j
            r = kron(r, x);
        else
            r = kron(r, y);
        end
    end
    if i == 1
        s = r;
    else
        s = s + r;
    end
end

end

