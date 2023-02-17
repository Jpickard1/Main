function xk = vecPower(x, k)
%VECPOWER  
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023
    xk = x;
    for i=1:k
        xk = kron(xk, x);
    end
end

