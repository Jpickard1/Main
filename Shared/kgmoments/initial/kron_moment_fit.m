function [rmat_params,objval,G] = kron_moment_fit(gdata,notris)
% KRON_MOMENT_FIT
%  [a,b,c,r] = kron_moment_fit(gdata)
%     gdata.nverts: the number of vertices
%     gdata.nedges: the number of undirected edges
%     gdata.nwedges: the number of wedges
%     gdata.ntripins: the number of tripins
%     gdata.ntris: the number of triangles

if ~exist('notris','var') || isempty(notris), notris=0; end

if ~isfield(gdata,'nverts')
    gdata.nverts = length(gdata.degs);
end

% Setup expectations
r = ceil(log2(gdata.nverts));
e = @(a,b,c,r) [
    (1/2)*((a+2*b+c)^r - (a+c)^r)
    (1/2)*(((a+b)^2 + (b+c)^2)^r - 2*(a*(a+b)+c*(c+b))^r - (a^2 + 2*b^2+c^2)^r + 2*(a^2 + c^2)^r)
    (1/6)*(((a+b)^3 + (b+c)^3)^r - 3*(a*(a+b)^2 + c*(b+c)^2)^r - 3*(a^3 + c^3 + b*(a^2+c^2)+b^2*(a+c) + 2*b^3)^r + 2*(a^3 + 2*b^3 + c^3)^r + 5*(a^3 + c^3 + b^2*(a+c))^r + 4*(a^3 + c^3 + b*(a^2 + c^2))^r - 6*(a^3 + c^3)^r)
    (1/6)*((a^3 + 3*b^2*(a+c) + c^3)^r - 3*(a*(a^2 + b^2) + c*(b^2 + c^2))^r + 2*(a^3+c^3)^r)];

f = zeros(4,1);
f(1) = gdata.nedges;
f(2) = gdata.nwedges;
f(3) = gdata.ntripins;
f(4) = gdata.ntris;

%% Minimize
%obj = @(x) sum( diag([1,1,1,1])*(f - e(x(1),x(2),x(3),r)).^2./(max(e(x(1),x(2),x(3),r),eps(1)) ));
%obj = @(x) sum( (diag([1,1,1,1])*((f - e(x(1),x(2),x(3),r))./(e(x(1),x(2),x(3),r))) ).^2 );
if notris
    weight=[1,1,1,0];
else
    weight=[1,1,1,1];
end
obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r))./(max(f,eps(1)) ));
confun = @(x) deal( x(3)-x(1) , [] ); % c-a <= 0 => c<=a
lb = [0 0 0];
ub = [1 1 1];
options = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Algorithm','active-set','Display','off');

minobj = Inf;
xmin = [];
for i=1:10
    x0 = sort(rand(1,3),2,'descend');
    x = fmincon(obj,x0,[],[],[],[],lb,ub,confun,options);
    %[x0; x]
    if obj(x)<minobj, xmin = x; minobj = obj(x); end
end    

%% Minimize (grid)
npts = 50;
[x1,x2,x3]=ndgrid(0:1/npts:1, 0:1/npts:1, 0:1/npts:1);
X=[x1(:), x2(:), x3(:)]';
x=gridmin(obj, X, 0);
if obj(x)<minobj, xmin = x; end



%% output
x = xmin;
a = x(1);
b = x(2);
c = x(3);
rmat_params=struct();
rmat_params.a = a;
rmat_params.b = b;
rmat_params.c = c;
rmat_params.r = r;
rmat_params.mat = [a b; b c];
objval = obj(x);

%% report
if nargout>1
    
    G = rmat(r,[a b; b c]); 
    G = triu(G,1); G=G|G';
    G = G - diag(diag(G));
    
    props = graph_properties(G);

    g = zeros(4,1); % f(1)=edges, f(2)=hairpins, f(3)=tripins, f(4)=triangles
    g(1) = props.nedges;
    g(2) = props.nwedges;
    g(3) = props.ntripins;
    g(4) = props.ntris;

    % report stats
    fprintf('\n');
    if isfield(gdata,'name')
        fprintf('graph               = %s\n', gdata.name);
    end
    fprintf('src vertices        = %20g\n', gdata.nverts);
    fprintf('objective value     = %20g\n', obj(x));
    fprintf('fitted params a,b,c = [%6.4f , %6.4f , %6.4f]\n', a,b,c);
    fprintf('\n');
    fprintf('%5s',''); fprintf(' %17s %17s %17s %17s','edges','hairpins','tripins','triangles'); fprintf('\n');
    fprintf('%5s','src'); fprintf(' %17i', f); fprintf('\n');
    fprintf('%5s','fit'); fprintf(' %17i', g); fprintf('\n');
    fprintf('%5s','est'); fprintf(' %17i', round(e(a,b,c,r))); fprintf('\n');
end
