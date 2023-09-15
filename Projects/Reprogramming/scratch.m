
clear

n = 4;
m = 3;
p = 10;

A = rand(n,m);
B = rand(n,p);
C = zeros(p,m);
S = randi([0,1],[p,m])

A/B

B*C

A(:,1) = B * C(:,1)

B \ A(:,1)
A(:,1) \ B

%%

clear

n = 300;
m = 40;
p = m;

B = rand(n,p);
Ctrue = full(sprand(p,m,0.2));
S = (Ctrue ~= 0);
A = B * Ctrue;

tic;
Clearned = structuredLS(A,B,S);
ts = toc;
tic;
Cdense = B\A;
td = toc;

nd = norm(Ctrue-Cdense);
ns = norm(Ctrue-Clearned);

disp("A: " + string(size(A)));
disp("B: " + string(size(B)));
disp("C: " + string(size(Ctrue)));
disp("Sparse Error: " + string(ns) + ", Time: " + string(ts));
disp("Dense Error: " + string(nd) + ", Time: " + string(td));


%% Can structuredLS compete with DMD?

clear;

n = 200;
t = 1000;

Atrue = full(sprand(n,n,0.05));
Atrue = Atrue ./ sum(Atrue,1);
% sum(Atrue,1)
X = zeros(n,t);
X(:,1) = rand(n,1);
for i=2:t
    X(:,i) = Atrue * X(:,i-1);
end

disp(find(isnan(X)));


% t = 20000;
% X = X(:,1:t);

tic;
out = DMD(X,[],1);
td = toc; disp('DMD');
tic;
As = SDMD(X,(Atrue ~= 0));
ts = toc; disp('Sparse');
tic;
Ae = exactDMD(X);
te = toc; disp('Exact');

nd = norm(Atrue - out.DMD.A_bar);
ns = norm(Atrue - As);
ne = norm(Atrue - Ae);

X1 = X(:,1:end-1);
X2 = X(:,2:end);

Xd2 = out.DMD.A_bar * X1;
ed = sum(abs(X2 - Xd2),'all');
Xs2 = As * X1;
es = sum(abs(X2 - Xs2),'all');
Xe2 = Ae * X1;
ee = sum(abs(X2 - Xe2),'all');

disp("X: " + string(size(X)));
disp("A: " + string(size(Atrue)));
disp("Sparse Error: " + string(ns) + ", Time: " + string(ts) + ", Prediction Error: " + string(es));
disp("DMD Error: " + string(nd) + ", Time: " + string(td) + ", Prediction Error: " + string(ed));
disp("Exact Error: " + string(ne) + ", Time: " + string(te) + ", Prediction Error: " + string(ee));

%%
St = (Atrue ~= 0)
Ss = (As ~= 0)
Sd = (out.DMD.A_bar > 1e-4)

Sd == St

AdmdS = zeros(size(out.DMD.A_bar));
AdmdS(out.DMD.A_bar > 1e-4) = out.DMD.A_bar(out.DMD.A_bar > 1e-4);
norm(Atrue - AdmdS)
