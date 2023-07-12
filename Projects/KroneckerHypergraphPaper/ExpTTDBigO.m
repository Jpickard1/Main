%% TTD Big-O Experiment
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

clear

tA = [];
tB = [];
tC = [];

maxN = 20^2;
for n=1:maxN
    if round(sqrt(n)) ~= sqrt(n)
        continue
    end
    disp(n)
    ns = sqrt(n);
    B = rand(ns,ns,ns);
    C = rand(ns,ns,ns);
    A = superkron(B,C);

    tic;
    eb = tt_tensor(B);
    tB = [tB; toc];
    tic;
    ec = tt_tensor(C);
    tC = [tC; toc];
    tic;
    ea = tt_tensor(A);
    tA = [tA; toc];
   
end

%% Make figure to describe the run time
figure; tiledlayout(1,2);

nexttile;
n = 31;
tt = tB+tC; tt = tt(1:n);
scatter([1:n],tA); hold on;
scatter([1:n],tt); hold on;
xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Time (sec.)','Interpreter','latex');
title('TTD Compute Time','Interpreter','latex');
legend(["Direct", "Kronecker"]);


nexttile;
X = [1:n];
d = 3;
r = 4;
D = 3 * (X .^ 2) * r^2;
K = 2 * 3 * X * r;

scatter(X, D); hold on;
scatter(X, K);
set(gca,'XTick',[])
set(gca,'YTick',[])
title('Theoretical Complexity', 'Interpreter', 'latex');
xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$', 'Interpreter', 'latex');
ylabel('Time', 'Interpreter', 'latex');
legend(["Direct", "Kronecker"]);


%% CPD calculations


clear; close all; clc;

tA = [];
tB = [];
tC = [];

errors = [];

maxN = 20^2;
r = 4;
for n=4:maxN
    if round(sqrt(n)) ~= sqrt(n)
        continue
    end
    disp(n)
    ns = sqrt(n);
    B = rand(ns,ns,ns);
    C = rand(ns,ns,ns);
    A = superkron(B,C);

    tic;
    Bc = cp_als(tensor(B), r);
    tB = [tB; toc];
    tic;
    Cc = cp_als(tensor(C), r);
    tC = [tC; toc];
    tic;
    Ac = cp_als(tensor(A), r^2);
    tA = [tA; toc];

    e1 = sum(abs(double(A - full(Ac))), 'all');
    Bc = double(full(Bc));
    Cc = double(full(Cc));
    e2 = sum(abs(A - superkron(Bc, Cc)), 'all');

    errors = [errors; e1 e2];
   
end

figure;
tiledlayout(1,2);
nexttile;
n = 19;
tt = tB+tC; tt = tt(1:n);
scatter([1:n],tA); hold on;
scatter([1:n],tt); hold on;
xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Time (sec.)','Interpreter','latex');
title('CPD Compute Time','Interpreter','latex');

nexttile;
plot(errors(:,1)); hold on;
plot(errors(:,2)); hold on;
xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Approximation Error','Interpreter','latex');
title('CPD Approximation Error','Interpreter','latex');
legend(["Direct", "Kronecker"]);
