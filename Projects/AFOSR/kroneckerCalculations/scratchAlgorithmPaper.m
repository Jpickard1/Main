%% Scratch Algorithm Paper
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 17, 2023

clear all; close all; clc;

n = 5;
A = rand(n^2,n^2,n^2);
[A1, A2] = NTKP(A);
tic; [B, Sigma] = tkpsvd(A, n * ones(1,6)); disp(toc);
tic; A12 = superkron(A1, A2); disp(toc);
B12 = Sigma(1) * superkron(B{1,1}, B{2,1});

disp(sum(sum(sum(abs(A - B12)))) / sum(sum(sum(abs(A)))))
disp(sum(sum(sum(abs(A - A12)))) / sum(sum(sum(abs(A)))))



%% Experiment

itrs = 5;
maxN = 20;
E1 = zeros(maxN, itrs);
T1 = zeros(maxN, itrs);
E2 = zeros(maxN, itrs);
T2 = zeros(maxN, itrs);

for n=2:maxN
    for i=1:itrs
        A = rand(n^2,n^2,n^2);
        [A1, A2] = NTKP(A);
        tic; [B, Sigma] = tkpsvd(A, n * ones(1,6), 1); 
        T1(n,i) = toc;
        tic; A12 = superkron(A1, A2); 
        T2(n,i) = toc;
        B12 = Sigma(1) * superkron(B{1,1}, B{2,1});
        
        E1(n,i) = (sum(sum(sum(abs(A - B12)))) / sum(sum(sum(abs(A)))));
        E2(n,i) = (sum(sum(sum(abs(A - A12)))) / sum(sum(sum(abs(A)))));
    end
    disp(n)
end

%% Plotting

N = (1:maxN) .^ 2;

figure;
subplot(1,2,1);
title('Error Comparison'); hold on;
xlabel('N'); ylabel('Error');
Y = mean(E2, 2);
STD = std(E2, 0, 2);
c1 = Y + STD;
c2 = Y - STD;
N2 = [N, fliplr(N)];
inBetween = [c1; fliplr(c2)];
fill(N2, inBetween, 'g');
plot(N, Y);
Y = mean(E1, 2);
STD = std(E1, 0, 2);
c1 = Y + STD;
c2 = Y - STD;
N2 = [N, fliplr(N)];
inBetween = [c1; fliplr(c2)];
fill(N2, inBetween, 'g');
plot(N, Y);
legend(["NTKP", "TKPSVD"]);
subplot(1,2,2,'YScale', 'log');
title('Time Comparison'); hold on;
xlabel('N'); ylabel('Time');
Y = mean(T2, 2);
STD = std(T2, 0, 2);
c1 = Y + STD;
c2 = Y - STD;
N2 = [N, fliplr(N)];
inBetween = [c1; fliplr(c2)];
fill(N2, inBetween, 'g');
plot(N, Y);
Y = mean(T1, 2);
STD = std(T1, 0, 2);
c1 = Y + STD;
c2 = Y - STD;
N2 = [N, fliplr(N)];
inBetween = [c1; fliplr(c2)];
fill(N2, inBetween, 'g');
plot(N, Y);
legend(["NTKP", "TKPSVD"]);

%% Old Plotting

figure; title('Error Comparison'); hold on;
xlabel('N'); ylabel('Error');
plot(mean(E2, 2));
plot(mean(E1, 2));
figure; title('Time Comparison'); hold on;
xlabel('N'); ylabel('Time');
plot(mean(T2, 2));
plot(mean(T1, 2));