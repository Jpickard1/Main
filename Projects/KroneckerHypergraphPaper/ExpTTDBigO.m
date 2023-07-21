%% TTD Big-O Experiment
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

% clear

n = 25;
itrs = 5;
tA = zeros(n, itrs);
tB = zeros(n, itrs);
tC = zeros(n, itrs);

for i=1:n
    for j=1:itrs
        disp(i);
        % ns = sqrt(i);
        ns = i;
        B = rand(ns,ns,ns);
        C = rand(ns,ns,ns);
        A = superkron(B,C);
    
        tic;
        eb = tt_tensor(B);
        tB(i,j) = toc;
        tic;
        ec = tt_tensor(C);
        tC(i,j) = toc;
        tic;
        ea = tt_tensor(A);
        tA(i,j) = toc;
    end
end

tt_TTD = tB+tC;
tA_TTD = tA;
%%

n = 500;
itrs = 3;

tA = zeros(n, itrs);
tB = zeros(n, itrs);
tC = zeros(n, itrs);
r = 4;

for i=1:n
    for j=1:itrs
        disp(i);
        ns = i;
        B = rand(ns,ns,ns);
        C = rand(ns,ns,ns);
        A = superkron(B,C);

        tic;
        Bc = cp_als(tensor(B), r);
        tB(i,j) = toc;
        tic;
        Cc = cp_als(tensor(C), r);
        tC(i,j) = toc;
        tic;
        Ac = cp_als(tensor(A), r^2);
        tA(i,j) = toc;
    end
end
tt_CPD = tB + tC;
tA_CPD = tA;

%% Make the finalized figure
set(groot, 'DefaultFigureRenderer', 'painters');
figure('Renderer', 'painters', 'Position', [0 0 900 400]);
tiledlayout(1,2);
fs = 18;

nexttile;
tt = tt_TTD; tA = tA_TTD;
% Calculate mean and standard error for each row
mean_tt = mean(tt, 2);
std_error_tt = std(tt, 0, 2) / sqrt(itrs);
mean_ta = mean(tA, 2);
std_error_ta = std(tA, 0, 2) / sqrt(itrs);

% Create scatter plots
% figure;
hold on;
errorbar(1:25, mean_ta, std_error_ta, 'ro', 'MarkerSize', 5, 'LineWidth', 1);
errorbar(1:25, mean_tt, std_error_tt, 'bo', 'MarkerSize', 5, 'LineWidth', 1);
hold off;

xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Time (sec.)','Interpreter','latex');
title('TTD Compute Time','Interpreter','latex');
% legend(["Direct", "Kronecker"]);
set(gca, 'YScale', 'log');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');

nexttile;
tt = tt_CPD; tA = tA_CPD;

% Calculate mean and standard error for each row
mean_tt = mean(tt, 2);
std_error_tt = std(tt, 0, 2) / sqrt(itrs);
mean_ta = mean(tA, 2);
std_error_ta = std(tA, 0, 2) / sqrt(itrs);

% Create scatter plots
% figure;
hold on;
errorbar(1:n, mean_ta, std_error_ta, 'ro', 'MarkerSize', 5, 'LineWidth', 1);
errorbar(1:n, mean_tt, std_error_tt, 'bo', 'MarkerSize', 5, 'LineWidth', 1);
hold off;

xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Time (sec.)','Interpreter','latex');
title('CPD Compute Time','Interpreter','latex');
legend(["Direct", "Kronecker"],'interpreter','latex');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');

% saveas(gcf, "tensorDecompositions_i5_v3.png")


%% CPD Experiment
% clear

n = 20;
itrs = 2;
tA = zeros(n, itrs);
tB = zeros(n, itrs);
tC = zeros(n, itrs);
r = 4;

for i=1:n
    for j=1:itrs
        disp(i);
        ns = i;
        B = rand(ns,ns,ns);
        C = rand(ns,ns,ns);
        A = superkron(B,C);

        tic;
        Bc = cp_als(tensor(B), r);
        tB(i,j) = toc;
        tic;
        Cc = cp_als(tensor(C), r);
        tC(i,j) = toc;
        tic;
        Ac = cp_als(tensor(A), r^2);
        tA(i,j) = toc;
    end
end

% Make the finalized figure
figure;
tt = tB+tC; tt = tt(1:n,:);

% Calculate mean and standard error for each row
mean_tt = mean(tt, 2);
std_error_tt = std(tt, 0, 2) / sqrt(itrs);
mean_ta = mean(tA, 2);
std_error_ta = std(tA, 0, 2) / sqrt(itrs);

% Create scatter plots
figure;
hold on;
errorbar(1:n, mean_ta, std_error_ta, 'ro', 'MarkerSize', 5, 'LineWidth', 1);
errorbar(1:n, mean_tt, std_error_tt, 'bo', 'MarkerSize', 5, 'LineWidth', 1);
hold off;

xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Time (sec.)','Interpreter','latex');
title('CPD Compute Time','Interpreter','latex');
legend(["Direct", "Kronecker"]);
set(gca, 'YScale', 'log');


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


% clear; close all; clc;

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
