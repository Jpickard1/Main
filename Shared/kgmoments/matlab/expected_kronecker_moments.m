function expdata = expected_kronecker_moments(params, r)
% EXPECTED_KRONECKER_MOMENTS
% 
% expdata = expected_kronecker_moments(params,r) takes a vector of
% parameters params = [a,b,c] and r, and returns the expected moment counts
% as follows:
%   expdata.nverts = 2^r
%   expdata.nedges = expected edges
%   expdata.nwedges = expected wedges
%   expdata.ntris = expected triangles
%   expdata.ntripins = expected tripins
%
% Alternatively, params can be a single struct output from
% kronecker_moment_fit, which has fields:
%   params.a, params.b, params.c, and params.r
% in which case, the value of r is unnecesssary, but it overrides the value
% in the parameter structure.
%
%

% David F. Gleich
% Sandia National Labs

% History
% :2011-04-12: Initial coding

if isstruct(params)
    a = params.a;
    b = params.b;
    c = params.c;
    
    if nargin==1
        r = params.r;
    else 
        % r is set by the command line argument.
    end
else
    assert(length(params) == 3);
    a = params(1);
    b = params(2);
    c = params(3);
    assert(nargin == 2); % make sure r got set
end

% formulas copied from the paper and kron


e = @(a,b,c,r) [
    (1/2)*((a+2*b+c)^r - (a+c)^r)
    (1/2)*(((a+b)^2 + (b+c)^2)^r - 2*(a*(a+b)+c*(c+b))^r - (a^2 + 2*b^2+c^2)^r + 2*(a^2 + c^2)^r)
    (1/6)*(((a+b)^3 + (b+c)^3)^r - 3*(a*(a+b)^2 + c*(b+c)^2)^r - 3*(a^3 + c^3 + b*(a^2+c^2)+b^2*(a+c) + 2*b^3)^r + 2*(a^3 + 2*b^3 + c^3)^r + 5*(a^3 + c^3 + b^2*(a+c))^r + 4*(a^3 + c^3 + b*(a^2 + c^2))^r - 6*(a^3 + c^3)^r)
    (1/6)*((a^3 + 3*b^2*(a+c) + c^3)^r - 3*(a*(a^2 + b^2) + c*(b^2 + c^2))^r + 2*(a^3+c^3)^r)];

moments = e(a,b,c,r);
expdata = [];
expdata.nverts = 2^r;
expdata.nedges = moments(1);
expdata.nwedges = moments(2);
expdata.ntripins = moments(3);
expdata.ntris = moments(4);
