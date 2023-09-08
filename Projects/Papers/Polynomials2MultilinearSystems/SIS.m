function [Am] = SIS(beta, gamma, stoch)
%SIS Constructs a matrix representation of the SIS system
%
% SIS SYSTEM
%  ds/dt = -beta S I + gamma I - stoch
%  di/dt =  beta S I - gamma I + stoch
%
% MATRIX REPRESENTATION
%     ss si is ii s i c
%      1  2  3  4 5 6 7
%
% REFERENCE: Kermack and McKendrick
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 17, 2023

if nargin == 2
    stoch = 0;
end

n = 2;
d = 2;
esv = numPolyTerms(n, d);
Am = sparse(n, esv);

Am(1,2) = - beta;   % ds/dt contains - beta S I
% Am(1,3) = - 0.5 * beta;   % ds/dt contains - beta S I
Am(1,6) =  gamma;         % ds/dt contains gamma I
Am(1,7) =  -stoch;        % ds/dt contains stoch
Am(2,2) =   beta;    % di/dt contains beta S I
% Am(2,3) =  0.5 * beta;    % di/dt contains beta S I
Am(2,6) =  -gamma;        % di/dt contains - gamma I
Am(2,7) =  stoch;         % di/dt contains - stoch

end

