function npt = numPolyTerms(nvars, degree)
%NUMPOLYTERMS Number of polynomial terms in the system
%   n - number of variables
%   d - maximum polynomial degree
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 6, 2023
    npt = 1;
    for i=1:degree
        npt = npt + nvars^i;
    end
end

