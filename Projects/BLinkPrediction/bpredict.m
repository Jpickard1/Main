function [Ei] = bpredict(V, E, n, b, f, Ec)
%BPREDICT Batched link prediction
%
%   This function performs batched link prediction according to the
%   algorithm in Joshua's research notebook from 03/09/2023
%
%   PARAMS:
%       V: vertex set
%       E: edge set
%       n: number of predictions to make
%       b: batch size
%       f: link prediction function
%       Ec: set of candidate edges
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 9, 2023

if nargin == 5
    % Set possible edges to predict
    Ec = nchoosek(1:V, 2);         % list all possible edges
    Ec = setdiff(Ec, E, 'rows');   % remove edges that already exist
end

Ek = E;

[~, R] = rankCandidateLinks(V,E,f,Ec);
i = 0;
while i < n
    % Update predicted edge set
    E = [E; Ec(R(1:b), :)];
    Ec(R(1:b), :) = [];

    % Update edge ranking
    [~, R] = rankCandidateLinks(V,E,f,Ec);

    % Update number of predicted new edges
    i = i + b;
end

Ei = setdiff(E, Ek, 'rows');

end
