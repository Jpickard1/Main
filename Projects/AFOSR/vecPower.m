function xk = vecPower(x, k)
%VECPOWER This function performs kronecker exponentiaion on a vector.
%
%PARAMETERS
%   x: vector
%   k: power
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023
    xk = x;
    for i=1:k-1
        xk = kron(xk, x);
    end
end

