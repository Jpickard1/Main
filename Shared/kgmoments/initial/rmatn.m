function A = rmatn(scale,P,noise,nedges)
% RMATN : generate power-law directed graph with R-MAT algorithm
%
% A = rmatn(scale)   returns a graph (adj matrix) with 2^scale vertices,
%                   with vertex numbers not randomized.
%
% A = rmatn(scale,P,noise)  return a graph (adj matrix) with 2^scale
% vertices, with vertex numbers not randomized.  The noise level controls
% what fraction of the maximum noise is added. A reasonable value is 0.5.


% Implementation of the Recursive MATrix (R-MAT) power-law graph
% generation algorithm (Chakrabati, Zhan & Faloutsos).
% This is a "nice" implementation in that  it is completely embarrassingly 
% parallel and does not require ever forming the adjacency matrix.
% Original by Jeremy Kepner, repackaged by John R. Gilbert, summer 2006
% Modified by David Gleich, 2009-02-03

% Changed to uniform random noise model following Seshadhri, Pinar, and
% Kolda (2011).

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
    P = [a b; c d];
    Ne = 8*Nv;
else
    Ne = ceil(((sum(sum(P)))^scale - (sum(diag(P)))^scale));
    P = P./sum(sum(P));
end

if exist('nedges','var')
    Ne = nedges*Nv;
end

maxnoise = min((P(1,1)+P(2,2))/2,P(1,2));
noiselvl = maxnoise*noise;

% Create index arrays.
ii = ones(Ne,1);
jj = ones(Ne,1);
% Loop over each order of bit.
ab = P(1,1)+P(1,2);
c_norm = P(2,1)./(P(2,1)+P(2,2));
a_norm = P(1,1)./(P(1,1)+P(1,2));
for ib = 1:lgNv
  % Compare with probabilities and set bits of indices.
  ii_bit = rand(Ne,1) > ab;
  jj_bit = rand(Ne,1) > ( c_norm.*ii_bit  + a_norm.*not(ii_bit) );
  ii = ii + (2^(ib-1)).*ii_bit;
  jj = jj + (2^(ib-1)).*jj_bit;
  %P = P + 1e-3*rand(size(P)); P = P./sum(sum(P));
  Pc = P;
  b = 2*noiselvl*rand(1); - noiselvl;
  Pc(1,2) = Pc(1,2) + b;
  Pc(2,1) = Pc(2,1) + b;
  f = (1-2*b)/(P(1,1) + P(2,2));
  Pc(1,1) = Pc(1,1)*f;
  Pc(2,2) = Pc(2,2)*f;
  a = Pc(1,1); b = Pc(1,2); c = Pc(2,1); d = Pc(2,2);
  ab = a+b;
  c_norm = c./(c+d);
  a_norm = a./(a+b);
end

A = sparse(ii,jj,1,Nv,Nv);
