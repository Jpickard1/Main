function rval=kron_moments(r, params, asarray)

a = params(1);
b = params(2);
c = params(3);

nedges = (1/2)*((a+2*b+c)^r - (a+c)^r);
nwedges = (1/2)*(((a+b)^2 + (b+c)^2)^r - 2*(a*(a+b)+c*(c+b))^r - (a^2 + 2*b^2+c^2)^r + 2*(a^2 + c^2)^r);
ntripins = (1/6)*(((a+b)^3 + (b+c)^3)^r - 3*(a*(a+b)^2 + c*(b+c)^2)^r - 3*(a^3 + c^3 + b*(a^2+c^2)+b^2*(a+c) + 2*b^3)^r + 2*(a^3 + 2*b^3 + c^3)^r + 5*(a^3 + c^3 + b^2*(a+c))^r + 4*(a^3 + c^3 + b*(a^2 + c^2))^r - 6*(a^3 + c^3)^r);
ntris = (1/6)*((a^3 + 3*b^2*(a+c) + c^3)^r - 3*(a*(a^2 + b^2) + c*(b^2 + c^2))^r + 2*(a^3+c^3)^r);


if nargin<=2 || asarray~=1
    rval.nverts = 2^r;
    rval.nedges = nedges;
    rval.nwedges = nwedges;
    rval.ntripins = ntripins;
    rval.ntris = ntris;
else 
    rval = [nedges; nwedges; ntripins; ntris]
end

