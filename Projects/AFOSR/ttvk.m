function pp = ttvk(T,v)
%TTVK Summary of this function goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023
    P = ttv(T,v,1);
    while ndims(P) > 1
        P = ttv(P,v,1);
    end
    P = tenmat(P, 1);
    pp = P(:);
end

