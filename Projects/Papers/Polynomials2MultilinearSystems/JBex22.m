function [Am] = JBex22()
%John Baillieul: Controllability and Observability of Polynomial Dynamical
% Systems Example 2.2
%
% EXAMPLE:
%   dx/dt = x^2
%   dy/dt = -x*y + u    // in this construction, u is added later in B
%
% MATRIX REPRESENTATION
%   xx xy yx yy x y c
%    1  2  3  4 5 6 7
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 1, 2023

n = 2;
d = 2;
esv = numPolyTerms(n, d);
Am = sparse(n, esv);
Am(1,4) = 1;
Am(2,2) = -1;

end
