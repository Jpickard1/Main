function [Am] = JCex2()
%JCex2 Jurdevic and Kupka Example 2
%
% EXAMPLE:
%   dx/dt = u           // in this construction, u is added later in B
%   dy/dt = x^2 + y^3
%
% MATRIX REPRESENTATION
%   xxx xxy xyx xyy yxx yxy yyx yyy xx xy yx yy  x  y  c
%     1   2   3   4   5   6   7   8  9 10 11 12 13 14 15
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 1, 2023

n = 2;
d = 3;
esv = numPolyTerms(n, d);
Am = sparse(n, esv);
Am(2,1) = 1;
Am(2,9) = 1;

end
