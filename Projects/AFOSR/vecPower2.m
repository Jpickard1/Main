function xk = vecPower2(x, s)
%VECPOWER This function performs kronecker exponentiaion on a vector until
%   the product contains s elements.
%
%PARAMETERS
%   x: vector
%   k: power
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 17, 2023
    xk = x;
    while numel(xk) < s
        xk = kron(xk, x);
    end
end

