function ob = obsvSym(a,c)
%OBSV  Compute the observability matrix.
%
%   OB = OBSV(A,C) returns the observability matrix [C; CA; CA^2 ...]
%
%   CO = OBSV(SYS) returns the observability matrix of the 
%   state-space model SYS with realization (A,B,C,D).  This is 
%   equivalent to OBSV(sys.a,sys.c).
%
%   For ND arrays of state-space models SYS, OB is an array with N+2
%   dimensions where OB(:,:,j1,...,jN) contains the observability 
%   matrix of the state-space model SYS(:,:,j1,...,jN).  
%
%   See also OBSVF, SS.

%   Thanks to Joseph C. Slater (Wright State University) and
%             Jesse A. Leitner (AFRL/VSSS, USAF)
%   Copyright 1986-2009 The MathWorks, Inc.
if nargin==1 && ~isnumeric(a)
   error(message('Control:general:NotSupportedModelsofClass','obsv',class(a)))
end
narginchk(2,2)

% Dimension checking
n = size(a,1);
ny = size(c,1);
if ~isequal(n,size(a,2),size(c,2)),
   error(message('Control:ltiobject:obsv1'))
end

% Allocate OB and compute each C A^k term
ob = sym('x',[n*ny,n]);
ob(1:ny,:) = c;
for k=1:n-1
  ob(k*ny+1:(k+1)*ny,:) = ob((k-1)*ny+1:k*ny,:) * a;
end


