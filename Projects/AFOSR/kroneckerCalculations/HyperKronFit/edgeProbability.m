function [p] = edgeProbability(n, theta, u, v)
% Evaluates equation 5.2
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

n0 = size(theta,1);
k = log(n) / log(n0);

p = 1;
for i=1:k
    i1 = mod(floor((u - 1)/ n0^(i-1)), n0) + 1;
    i2 = mod(floor((v - 1)/ n0^(i-1)), n0) + 1;
    p = p * theta(i1, i2);
end


end