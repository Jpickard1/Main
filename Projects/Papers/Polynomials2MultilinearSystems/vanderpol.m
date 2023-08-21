function [Am] = vanderpol(mu)
%VANDERPOL Constructs a matrix representation of the Van Der Pol system
%
% CLASSIC VAN DER POL
%
%   (d2x/dt2) - mu (1-x^2) (dx/dt) + x = 0
%
% VAN DER POL EXPANSION
%
%   dx/dt = mu ( x - (1/3)x^3 - y)
%   dy/dt = (1/mu)x
%
% MATRIX REPRESENTATION
%   xxx xxy xyx xyy yxx yxy yyx yyy xx xy yx yy  x  y  c
%     1   2   3   4   5   6   7   8  9 10 11 12 13 14 15
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 7, 2023

n = 2;
d = 3;
esv = numPolyTerms(n, d);
Am = sparse(n, esv);

Am(1,13) = mu;
Am(1,1)  = - mu / 3;
Am(1,14) = - mu;
Am(2,13) = 1 / mu;

end

