function [Am] = lorenz(sigma, rho, beta)
%LORENZ Constructs a matrix representation of the Lorenz system as a
% polynomial
%
% LORENZ SYSTEM
%     dx/dt = sigma*y - sigma * x
%     dy/dt = -x*z +x*rho -y
%     dz/dt = x*y - beta * z
%
% MATRIX REPRESENTATION
%     xx xy xz yx yy yz zx zy zz  x  y  z  c
%      1  2  3  4  5  6  7  8  9 10 11 12 13
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 7, 2023

n = 3;
d = 2;
esv = numPolyTerms(n, d);
Am = sparse(n, esv);

Am(1,10) = -sigma;  % dx/dt contains -sigma*x
Am(1,11) =  sigma;  % dx/dt contains  sigma*x
Am(2,3)  = -1;      % dy/dt contains -x*z
Am(2,10) =  rho;    % dy/dt contains x*rho
Am(2,11) = -1;      % dy/dt contains -y
Am(3,2)  =  1;      % dz/dt contains x*y
Am(3,12) = -beta;   % dz/dt contains -beta*z

end

