function [EE] = allDirectedHyperedges(varargin)
%UNDIRECTEDHYPEREDGES
%
%   E is a m x k matrix for a k-uniform hypergraph on m vertices. Each
%   hyperedge is given in a particular order. The purpose of this function
%   is to return all possible orders of vertices in each hyperedge.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 15, 2023

if nargin == 1
    E = varargin{1};    % Get parameter
elseif nargin == 2
    n = varargin{1};    % Get parameter
    k = varargin{2};    % Get parameter
    A = zeros(n * ones(1,k));
    idxs = cell(1,k);
    [idxs{:}] = find(A==0);
    E = cell2mat(idxs);
else
    error(['Joshua, you wrote a bad function. Please fix this and stay' ...
        'compatable with allUndirectedHyperedgs.m']);
end

E = unique(sort(E, 2), 'rows');
k = size(E, 2);
rpts = factorial(k);
EE = zeros(size(E,1) * rpts, k);

for e=1:size(E,1)
    hedge = E(e,:);
    EE(((e-1) * rpts + 1):(e * rpts),:) = perms(hedge);
end

EE = unique(EE, 'rows');


end