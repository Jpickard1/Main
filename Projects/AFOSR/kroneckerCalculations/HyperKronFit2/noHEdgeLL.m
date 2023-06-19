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
elseif nargin == 3
    n = varargin{1};
    theta = varargin{2};
    idxs = varargin{3};
    LL = hedgeLL(n, theta, idxs);
else
    error('noHEdgeLL.m: the passed arguments are not accepted.')
end

nLL = log(1 - exp(LL));

end