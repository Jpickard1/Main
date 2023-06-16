function [LL] = hedge3LLapx(n, theta, u, v, w)
%HEDGE3LLAPX This function is one optimization of hedgeLLapx for 3-uniform
% hypergraphs.
%
% PARAMETERS:
%   n: number of vertices in Kronecker hypergraph
%   theta: Kronecker initiator matrix
%   u,,v,w: indices into Kronecker hypergraph
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 16, 2023

n0 = size(theta,1);
kronExp = log(n) / log(n0);

n0Exp = n0 .^ (1:kronExp-1);

LL = 0;
for i=1:kronExp
    i1 = mod(floor((u - 1) / n0Exp(i)), n0) + 1;
    i2 = mod(floor((v - 1) / n0Exp(i)), n0) + 1;
    i3 = mod(floor((w - 1) / n0Exp(i)), n0) + 1;
    LL = LL - theta(i1, i2, i3) - 0.5 * (theta(i1, i2, i3)^2);
end


end