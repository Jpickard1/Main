function [LL] = edgeLL(n, theta, u, v)
% Evaluates equation 5.2
%
%   u is row index
%   v is column index
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 6, 2023

% theta = theta / sum(sum(theta));

n0 = size(theta,1);
kronExp = log(n) / log(n0);

LL = 0;
for i=1:kronExp
    i1 = mod(floor((u - 1)/ n0^(i-1)), n0) + 1;
    i2 = mod(floor((v - 1)/ n0^(i-1)), n0) + 1;
    LL = LL + log(theta(i1, i2));
end


end