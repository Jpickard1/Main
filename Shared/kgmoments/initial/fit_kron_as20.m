%% Fit a stochastic Kronecker model

addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl-4.0');
addpath('~/dev/libbvg');

%%
A = readSMAT('~/data/leskovec/as20.smat');
A = A|A';
A = largest_component(A);
%%
% try with and without this section
A = A-diag(diag(A));
d = full(sum(A,2)-diag(A));
%% Compute features
d = full(sum(A,2));
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
obj = @(x) sum( diag([1,1,1,1])*(f - e(x(1),x(2),x(3),r)).^2./max(e(x(1),x(2),x(3),r),eps(1)) );
%obj = @(x) sum( (diag([1,1,1,1])*((f - e(x(1),x(2),x(3),r))./(e(x(1),x(2),x(3),r))) ).^2 );
%obj = @(x) sum( diag((sum(f)-f)/sum(f))*(f - e(x(1),x(2),x(3),r)).^2/max(sum(e(x(1),x(2),x(3),r)),eps(1)) );
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
%%
x= [1         0.722037778924319                         0];
%%
x= [1         0.721904855336955                         0];

%% Minimize (grid)
r = ceil(log2(num_vertices(A)));
obj = @(x) sum( diag([1,1,1,1])*(f - e(x(1),x(2),x(3),r)).^2./max(e(x(1),x(2),x(3),r),eps(1)) );
%obj = @(x) sum( diag((sum(f)-f)/sum(f))*(f - e(x(1),x(2),x(3),r)).^2/max(sum(e(x(1),x(2),x(3),r)),eps(1)) );
npts = 200;
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
% Jure's values
x = [0.992 0.57 0.05];

%% 
% Jure's values from Art's method
x = [0.99 0.73 0.01];

%% 
% Optimum over 
x = [0.99 0.73 0.00];

%%
a = x(1);
b = x(2);
c = x(3);
[f e(a,b,c,r)]
obj(x)
G = rmatn(r,[a b; b c]); 
G = triu(G,1); G=G|G';
G = G - diag(diag(G));
G = G(sum(G,2)>0,sum(G,2)>0);
dg = full(sum(G,2));

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
[N,X]=hist(dg,max(dg)); loglog(X,N,'r.');
[N,X]=hist(d,max(d)); loglog(X,N,'b.');
set(gca,'YScale','log'); set(gca,'XScale','log'); hold off;
title('degree');