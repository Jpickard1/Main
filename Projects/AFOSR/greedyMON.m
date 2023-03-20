function [D, OD] = greedyMON(O, n)
%GREEDYMON This function performs greedy node selection to get a set of
%   minimal observable nodes.
%
%   PARAMETERS:
%       - O: Observability matrices of individual vertices
%       - n: number of vertices
%   RETURNS:
%       - D: set of oberverd vertices
%       - OD: observability matrix
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 3, 2023

% Greedy Node Selection for Minimum Observable Nodes
S = 1:n;                                    % Unobserved vertices
D = [];                                     % Observed vertices
OD = [];                                    % Observability Matrix
while rank(OD) < n
    deltaS = zeros(length(S),1);              % Vector to store all changes in rank
    for i=1:length(S)                       % Try all vertices in S
        vx = S(i);                          % Get vertex
        ODS = [OD; O{vx}];                  % Set possible new observabliity matrix
        deltaS(i) = rank(ODS) - rank(OD);      % Compute improved rank for vx
    end
    [~, vx] = max(deltaS);                  % Greedy selection of vertex
    D = [D S(vx)];                          % Add new vertex to observe
    OD = [OD; O{S(vx)}];                    % Set new observability matrix
    S = S(S ~= vx);                         % Remove vx from S
end

end

