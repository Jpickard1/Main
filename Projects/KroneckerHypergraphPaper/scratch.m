clear all; close all; clc

n = 2;
itrs = 100;
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
xlabel("Time");
ylabel("Frequency");
title("Run Time");
