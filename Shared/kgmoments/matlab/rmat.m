function A = rmat(P,scale,nedges)
% RMAT : generate power-law directed graph with R-MAT algorithm
%
% A = rmat(T,scale)   returns a graph (adj matrix) with 2^scale vertices,
%                     with vertex numbers not randomized.


% Implementation of the Recursive MATrix (R-MAT) power-law graph
% generation algorithm (Chakrabati, Zhan & Faloutsos).
% This is a "nice" implementation in that  it is completely embarrassingly 
% parallel and does not require ever forming the adjacency matrix.
% Original by Jeremy Kepner, repackaged by John R. Gilbert, summer 2006
% Modified by David Gleich, 2009-02-03

% Set number of vertices.
lgNv = scale;
Nv = 2^lgNv;

% Set R-MAT probabilities.
% Create a single parameter family.
% Parameters can't be symmetric in order
% for it to be a power law distribution.
if ~exist('P','var') || isempty(P),
    p = 0.6; 
    a = p;  b = (1 - a)/3;  c = b;  d = b;
    Ne = 8*Nv;
else
    Ne = ceil(((sum(sum(P)))^scale - (sum(diag(P)))^scale));
    Ne = Ne + round(sqrt(Ne)*randn(1)); % throw a little bit of variance in the edges
    P = P./sum(sum(P));
    a = P(1,1); b = P(1,2); c = P(2,1); d = P(2,2);
end

if exist('nedges','var')
    Ne = nedges*Nv;
end

fprintf('Generating Kronecker model with %i vertices and %i edges\n', ...
    Nv, Ne);


[ii,jj] = rmatedges(a,b,c,d,Ne,lgNv);

A = sparse(ii,jj,1,Nv,Nv);

if nnz(A) < Ne
    [ii2,jj2] = rmatedges(a,b,c,d,Ne-nnz(A),lgNv);
    A = spones(A + sparse(ii2,jj2,1,Nv,Nv));
end


function [ii,jj] = rmatedges(a,b,c,d,Ne,lgNv)
% Create index arrays.
ii = ones(Ne,1);
jj = ones(Ne,1);
% Loop over each order of bit.
ab = a+b;
c_norm = c./(c+d);
a_norm = a./(a+b);
for ib = 1:lgNv
  % Compare with probabilities and set bits of indices.
  ii_bit = rand(Ne,1) > ab;
  jj_bit = rand(Ne,1) > ( c_norm.*ii_bit  + a_norm.*not(ii_bit) );
  ii = ii + (2^(ib-1)).*ii_bit;
  jj = jj + (2^(ib-1)).*jj_bit;
end


