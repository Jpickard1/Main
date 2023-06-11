%% Fit a stochastic Kronecker model

addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl-4.0');
addpath('~/dev/libbvg');

%% Load a graph

%% 
% Test graph
A = sparse(ones(32));

%%
G = bvgraph('/home/dgleich/data/webgraphs/itgraphs/wb-cs.stanford');
A = sparse(G);
A = A|A';
A = largest_component(A);
name ='wb-cs-stanford';
%%
A = readSMAT('/home/dgleich/data/router/itdk0304.smat');
A = A|A';
A = largest_component(A);
name ='itdk0304';
%%
A = readSMAT('~/data/arxiv/hep-th.smat');
A = A|A';
A = largest_component(A);
name ='arxiv-hep-th';
%%
A = readSMAT('~/data/leskovec/as20.smat');
A = A|A';
A = largest_component(A);
name ='as20';
%%
G = bvgraph('/home/dgleich/data/wikipedia/20051105/enwiki-20051105');
A = sparse(G);
A = A|A';
A = largest_component(A);
name ='enwiki-2005';

%% Normalize graph
A = A-diag(diag(A));
d = full(sum(A,2)-diag(A));
%% Compute features
f = zeros(4,1); % f(1)=edges, f(2)=hairpins, f(3)=tripins, f(4)=triangles
f(1) = sum(d)/2;
f(2) = sum(d.*(d-1))/2;
f(3) = sum(d.*(d-1).*(d-2))/6;
ccfs = clustering_coefficients(A,'undirected',1','unweighted',1);
tris = ccfs.*d.*(d-1);
f(4) = sum(tris)/6;

%% Setup expectations
e = @(a,b,c,r) [
    (1/2)*((a+2*b+c)^r - (a+c)^r)
    (1/2)*(((a+b)^2 + (b+c)^2)^r - 2*(a*(a+b)+c*(c+b))^r - (a^2 + 2*b^2+c^2)^r + 2*(a^2 + c^2)^r)
    (1/6)*(((a+b)^3 + (b+c)^3)^r - 3*(a*(a+b)^2 + c*(b+c)^2)^r - 3*(a^3 + c^3 + b*(a^2+c^2)+b^2*(a+c) + 2*b^3)^r + 2*(a^3 + 2*b^3 + c^3)^r + 5*(a^3 + c^3 + b^2*(a+c))^r + 4*(a^3 + c^3 + b*(a^2 + c^2))^r - 6*(a^3 + c^3)^r)
    (1/6)*((a^3 + 3*b^2*(a+c) + c^3)^r - 3*(a*(a^2 + b^2) + c*(b^2 + c^2))^r + 2*(a^3+c^3)^r)];

%% Minimize
r = ceil(log2(num_vertices(A)));
obj = @(x) sum( diag([1,1,1,1])*(f - e(x(1),x(2),x(3),r)).^2./(max(e(x(1),x(2),x(3),r),eps(1)) ));
%obj = @(x) sum( (diag([1,1,1,1])*((f - e(x(1),x(2),x(3),r))./(e(x(1),x(2),x(3),r))) ).^2 );
confun = @(x) deal( x(3)-x(1) , [] ); % c-a <= 0 => c<=a
lb = [0 0 0];
ub = [1 1 1];
options = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Algorithm','active-set');

for i=1:10
    x0 = sort(rand(1,3),2,'descend');
    x = fmincon(obj,x0,[],[],[],[],lb,ub,confun,options);
    [x0; x]
end    
npts = 100;

%% Minimize (grid)
r = ceil(log2(num_vertices(A)));
obj = @(x) sum( (f -e(x(1),x(2),x(3),r)).^2./max(e(x(1),x(2),x(3),r),eps(1)) );
npts = 100;
[x1,x2,x3]=ndgrid(0:1/npts:1, 0:1/npts:1, 0:1/npts:1);
X=[x1(:), x2(:), x3(:)]';
x=gridmin(obj, X);

%% 
% Pick the 3rd element over a finer range
obj3 = @(x3) obj([x(1) x(2) x3]);
x3min = gridmin(obj3, 0:1/(npts^2-1):1);
x3min

%% 
% Pick the 2n element over a finer range
obj2 = @(x2) obj([x(1) x2 x(3)]);
x2min = gridmin(obj2, 0:1/(npts^2-1):1);
x2min

%% 
% Pick the 1st element over a finer range
obj1 = @(x1) obj([x1 x(2) x(3)]);
x1min = gridmin(obj1, 0:1/(npts^2-1):1);
x1min

%%
% Assign from x3min
x(3) = x3min;

%%
% Assign from x2min
x(2) = x2min;

%%
% Assign from x1min
x(1) = x1min;

%% Construct the new graph
a = x(1);
b = x(2);
c = x(3);
G = rmat(r,[a b; b c]); 
G = triu(G,1); G=G|G';
G = G - diag(diag(G));
dg = full(sum(G,2)-diag(G));

g = zeros(4,1); % f(1)=edges, f(2)=hairpins, f(3)=tripins, f(4)=triangles
g(1) = sum(dg)/2;
g(2) = sum(dg.*(dg-1))/2;
g(3) = sum(dg.*(dg-1).*(dg-2))/6;
ccfsg = clustering_coefficients(G,'undirected',1','unweighted',1);
trisg = ccfsg.*dg.*(dg-1);
g(4) = sum(trisg)/6;

%%
% report stats
fprintf('\n');
fprintf('graph               = %s\n', name);
fprintf('src vertices        = %20g\n', size(A,1));
fprintf('objective value     = %20g\n', obj(x));
fprintf('fitted params a,b,c = [%6.4f , %6.4f , %6.4f]\n', a,b,c);
fprintf('\n');
fprintf('%5s',''); fprintf(' %17s %17s %17s %17s','edges','hairpins','tripins','triangles'); fprintf('\n');
fprintf('%5s','src'); fprintf(' %17i', f); fprintf('\n');
fprintf('%5s','fit'); fprintf(' %17i', g); fprintf('\n');
fprintf('%5s','est'); fprintf(' %17i', round(e(a,b,c,r))); fprintf('\n');
%%

figure(1);clf; hold on; loglog(-sort(-dg),'r'); loglog(-sort(-d),'b'); 
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('degree');

%%

figure(2); clf; hold on;
loglog(-sort(-(dg.*(dg-1))),'r'); loglog(-sort(-(d.*(d-1))),'b'); 
set(gca,'YScale','log'); set(gca,'XScale','log'); 
hold off; title('hairpins');

%%
figure(3); clf; hold on;
loglog(-sort(-(dg.*(dg-1).*(dg-2))),'r'); loglog(-sort(-(d.*(d-1).*(d-2))),'b'); 
set(gca,'YScale','log'); set(gca,'XScale','log'); 
hold off; title('tripins');

%%
trisg = clustering_coefficients(G,'undirected',0,'unweighted',0).*dg.*(dg-1);
figure(4);clf; hold on; loglog(dg,trisg,'r.'); loglog(d,tris,'b.'); 
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('triangles'); drawnow

%%
figure(1); clf; hold on;
[N,X]=hist(dg,min(max(dg),1e6)); loglog(X,N,'r.');
[N,X]=hist(d,min(max(d),1e6)); loglog(X,N,'b.');
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('degree'); xlabel('Degree at Vertex');  ylabel('Count');
print('-dpng',sprintf('%s-deg.png',name));

%%
figure(2); clf; hold on;
[N,X]=hist(dg.*(dg-1),min(max(dg.*(dg-1)),1e6)); loglog(X,N,'r.');
[N,X]=hist(d.*(d-1),min(max(d.*(d-1)),1e6)); loglog(X,N,'b.');
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('hairpins'); xlabel('Hairpins at Vertex');  ylabel('Count');

%%
figure(3); clf; hold on;
[N,X]=hist(dg.*(dg-1).*(dg-2),min(max(dg.*(dg-1).*(dg-2)),1e6)); loglog(X,N,'r.');
[N,X]=hist(d.*(d-1).*(d-2),min(max(d.*(d-1).*(d-2)),1e6)); loglog(X,N,'b.');
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('tripins'); xlabel('Tripins at Vertex'); ylabel('Count');
print('-dpng',sprintf('%s-tripins.png',name));

%%
figure(4); clf; hold on;
 loglog(d,ccfs,'b.'); loglog(dg,ccfsg,'r.');
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('triangles'); xlabel('Degree'); ylabel('Clustering Coefficient');
print('-dpng',sprintf('%s-ccfs.png',name));
