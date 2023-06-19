function [nLL] = noHEdgeLL(varargin)
%NOHEDGELL evaluates the log likelihood of an edge not existing in a hypergraph
%
% ARGUMENTS
%   n is size of Kronecker hypergraph
%   theta is Kronecker inidiator
%   idxs is the node in the Kronecker graph
%
% TODO
%   this is only configured for graphs. See hedgeLL.m
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

if nargin == 1
    LL = varargin{1};
if nargin == 3
    theta = varargin{1};
    idxs = varargin{2};
    LL = hedgeLL(theta, idxs);
end

nLL = log(1 - exp(LL));

end