function A = rmat2(P,r,nedges,noise)
% RMAT : generate power-law directed graph with R-MAT algorithm
%
% A = rmat(P,r)   returns a graph (adj matrix) with 2^r vertices,
%                 with vertex numbers not randomized.
% Implementation of the Recursive MATrix (R-MAT) power-law graph
% generation algorithm (Chakrabati, Zhan & Faloutsos).

% Original by Jeremy Kepner, repackaged by John R. Gilbert, summer 2006
% Modified by David Gleich, 2009-02-03, and again on 2009-02-13 to
% implement the algorithm from the paper

% Set R-MAT probabilities.
% Create a single parameter family.
% Parameters can't be symmetric in order
% for it to be a power law distribution.
if ~exist('nedges','var') || isempty(nedges)
    nedges = ceil(((sum(sum(P)))^r - (sum(diag(P)))^r));
end

if ~exist('noise','var') || isempty(noise), noise=1e-3; end

A = rmat2_impl(P,r,nedges,noise);
