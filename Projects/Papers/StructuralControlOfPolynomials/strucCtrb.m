function [ctrb] = strucCtrb(HG, CTRLS)
%STRUCCTRB Verifies if the polynomial defined by the directed hypergraph is
% controllable.
%
%   INPUTS:
%   * HG is a hypergraph struct from readEdgeSet
%   * CTRLS is a vector of controllable vertices
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

n = numel(HG.V);

%% Verify if all vertices are accessible from the CTRLS
A = HG2Aten(HG);
A = tensor(A);
X0 = zeros(n,1);
X0(CTRLS) = 1;

% CTRLS = [1 2]

VISITED = zeros(n,1);
VISITED(CTRLS) = 1;
i = 1;
while i < n && sum(VISITED) < n
    X0 = ttvk2(A,X0);
    VISITED(X0 ~= 0) = 1;
    i = i + 1;
end

if sum(VISITED) < n
    disp('Not Controllable. Not all vertices are accessible');
end

%% Check for dilations


end

