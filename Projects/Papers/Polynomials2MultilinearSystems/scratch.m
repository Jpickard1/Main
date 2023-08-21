%% Polynomials to Multilinear Systems
%
%   POLYNOMIAL REPRESENTATION: Amat * (x kron ... kron x)
%
%   MULTILINEAR REPRESENTATION: sptensor class by Kolda and Bader in tensor
%   toolbox
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 5, 2023

%% mvpoly construction
clear; close all; clc;

Am = rand(3,13);
p = mvpoly("Am",Am,"maxD",2);
x = rand(3,1);
sv = [kron(x,x);x;1];

disp(Am * sv)
disp(p.eval(x))

%% multilinear system construction from a polynomial
clear; close all; clc;

n = 2;
Am = rand(2,3);
p = mvpoly("Am",Am,"maxD",1);
m = multilinearSystem("poly",p);

x = rand(n,1);

p.eval(x)
m.eval(x)

%% Test Evaluation of multlinear system vs a polynomial
clear; close all; clc;

n = randi([2 5], 1);
d = randi([2 5], 1);
npt = numPolyTerms(n,d);
Am = rand(n, npt);
p = mvpoly("Am",Am,"maxD",d);
m = multilinearSystem("poly",p);

x = rand(n,1);

yp = p.eval(x);
ym = m.eval(x);

assert(sum(abs(yp - ym)) < 1e-10)

%% Test Evaluation of multlinear system vs a polynomial with sparse storage of polynomial
clear; close all; clc;

n = randi([2 10], 1);
d = randi([2 10], 1);
n = 5;
d = 10;
npt = numPolyTerms(n,d);
Am = sprand(n, npt, 0.01);
p = mvpoly("Am",Am,"maxD",d);
m = multilinearSystem("poly",p);

x = rand(n,1);

tic;
yp = p.eval(x);
tp = toc
tic;
ym = m.eval(x);
tm = toc

assert(sum(abs(yp - ym)) < 1e-10)

%% Lorenz System
n = 3; % x=1, y=2, z=3
d = 2;
esv = numPolyTerms(n, d);

sigma = 10;
beta = 8/3;
rho = 28;

Am = lorenz(sigma, rho, beta);

XX = sym('x',[3,1]);
disp(Am * [kron(XX,XX); XX; 1]);

p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);

[~,a] = ode45(@(t, x) p.eval(x),[0 100],[1 1 1]);
figure; plot3(a(:,1),a(:,2),a(:,3))

m = multilinearSystem('poly',p);

[~,a] = ode45(@(t, x) m.eval(x),[0 100],[1 1 1]);
figure; plot3(a(:,1),a(:,2),a(:,3))

%% Van Der Pol System

p = mvpoly('type','van der pol','mu',0.9);
[~,a] = ode45(@(t, x) p.eval(x),[0 1000],rand(2,1));
figure; plot(a(:,1),a(:,2)); title(p.title(), 'Interpreter','latex');



%%

sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
[t,a] = ode45(f,[0 100],[1 1 1]);     % Runge-Kutta 4th/5th order ODE solver
figure; plot3(a(:,1),a(:,2),a(:,3))


%%

for i=1:10000
    yr = rand(3,1);
    assert(sum(abs(f(0, yr) - p.eval(yr))) < 1e-5)
end

%%

delta = 1;
T = 1000;
X = zeros(n,T);
X(:,1) = rand(n,1);
for t=2:T
    X(:,t) = X(:,t-1) + delta * p.eval(X(:,t-1));
end

figure;
plot3(X(:,1), X(:,2), X(:,3))


%%
X = sptenrand([10 100],0.01)

reshape(X, [10 10 10])



%% Kronecker indices

% Example data
x = [1; 2; 3];
d = 4;
i = 7;
n = 3;
% Call the function to get the indices j_1, j_2, ..., j_d
indices = kroneckerIndices(n, d, 80);

disp(indices)

x = sym('x_%d',[3,1]);
KroneckerPower(x,4)

%% Subtensors of sparse tensors

% Create the 3x3x3x3 sparse tensor A
sizeA = [3, 3, 3, 3];
A = sptensor(sizeA);

% Create the 3x3x3 sparse tensor B
sizeB = [3, 3, 3];
B = sptensor(sizeB);
B(1,1,1) = 1;
B(2,2,2) = 2;
B(3,3,3) = 3;

% Set B as a sub-tensor of A at the first dimension (A(:,:,:,1) = B)
A(:,:,:,1) = B;


