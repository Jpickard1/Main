%% tensor eigenvalue figure

clear; close all; clc;
figure; tiledlayout(1,2);

nexttile; load("eigencalculationRunTimes.mat");
histogram(tA); hold on;
histogram(tB + tC);
xlabel("Time (sec.)");
ylabel("Frequency");
title("H-Eigenvalues of $A\in\mathbf{R}^{4\times 4\times 4}$",'interpreter','latex');
legend(["Standard", "Kronecker"])

nexttile; load("eigencalculationRunTimes6.mat");
histogram(tA); hold on;
histogram(tB + tC);
xlabel("Time (sec.)");
ylabel("Frequency");
title("H-Eigenvalues of $A\in\mathbf{R}^{6\times 6\times 6}$",'interpreter','latex');
legend(["Standard", "Kronecker"])



%% tensor eigenvalues
clear all; close all; clc

n = 2;
itrs = 1000;
tB = zeros(itrs,1);
tC = zeros(itrs,1);
tA = zeros(itrs,1);
BCAerror = zeros(itrs,1);
for i=1:itrs
    disp(i)
    B = rand(n,n,n);
    C = rand(n,n,n);
    A = superkron(B,C);
    
    tic;
    eb = heig(B);
    tB(i) = toc;
    tic;
    ec = heig(C);
    tC(i) = toc;
    tic;
    ea = heig(A);
    tA(i) = toc;

    BCAerrors(i) = max(ea) - max(eb)*max(ec);

end

figure;
histogram(abs(BCAerrors));
xlabel('Error');
ylabel('Frequency');
title('Calculation Error');

figure;
histogram(tA); hold on;
histogram(tB + tC);
xlabel("Time (sec.)",'Interpreter','latex');
ylabel("Frequency",'Interpreter','latex');
title("Time to Calculate Largest H-Eigenvalue of $$\textsf{A}\in\mathbf{R}^{4\times4\times4}$$",'Interpreter','latex');
legend(["Direct", "Kronecker"])
xlim([1.25, 1.8])

% save("eigencalculationRunTimes.mat")

max(max(BCAerrors))

%% 
%% tensor eigenvalues
% clear all; close all; clc

n1 = 2;
n2 = 3;
itrs = 100;
tB = zeros(itrs,1);
tC = zeros(itrs,1);
tA = zeros(itrs,1);
BCAerror = zeros(itrs,1);
for i=1:itrs
    disp(i)
    B = rand(n1,n1,n1);
    C = rand(n2,n2,n2);
    A = superkron(B,C);
    
    tic;
    eb = heig(B);
    tB(i) = toc;
    tic;
    ec = heig(C);
    tC(i) = toc;
    tic;
    ea = heig(A);
    tA(i) = toc;

    BCAerrors(i) = max(ea) - max(eb)*max(ec);

end

% save("eigencalculationRunTimes6.mat")

%% tensor contractions
clear all; close all; clc

n = 20;
itrs = 1000;
tB = zeros(itrs,1);
tC = zeros(itrs,1);
tA = zeros(itrs,1);
BCAerror = zeros(itrs,1);
for i=1:itrs
    disp(i)
    B = rand(n,n,n);
    C = rand(n,n,n);
    A = superkron(B,C);

    x = rand(n,1);
    y = rand(n,1);
    xy = kron(x,y);
    
    tic;
    vb = tensorprod(B, x, 1, 1);
    tB(i) = toc;
    tic;
    vc = tensorprod(C, y, 1, 1);
    vbc = kron(vb, vc);
    tC(i) = toc;
    tic;
    va = tensorprod(A, xy, 1, 1);
    tA(i) = toc;

    BCAerrors(i) = norm(va - vbc);

end

figure;
histogram(abs(BCAerrors));
xlabel('Error');
ylabel('Frequency');
title('Calculation Error');

figure;
histogram(tA, 10); hold on;
histogram(tB + tC, 10);
xlabel("Time (sec.)");
ylabel("Frequency");
title("Run Time");
legend(["Base Calculation", "Calculation with Kronecker"])

% save("contractionCalculationRunTimes.mat")

max(max(BCAerrors))