function pp = ttvk2(T,v, ndim)
%TTVK Summary of this function goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023
    P = ttv(T,v,1);
    while ndims(P) > ndim
        P = ttv(P,v,1);
    end
    P = tenmat(P, 1);
    pp = P(:);
end

