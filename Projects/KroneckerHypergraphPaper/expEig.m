%% H-eigenvalue calculations
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 1, 2023

% addpath(genpath('C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\AFOSR\'))
% addpath(genpath('C:\Users\picka\Documents\my_projects\DBTM\teneig-2.0 (1)\teneig-2.0-JP\'))

maxN = 35;
itrs = 5;

tA = zeros(maxN, itrs);
tBC = zeros(maxN, itrs);


for n=1:maxN
    disp(n);
    for i=1:itrs
        B = randSym3way(n);
        C = randSym3way(n);
        A = superkron(B,C);
        
        xB = rand(n, 1);
        xC = rand(n, 1);
        xA = kron(xB,xC);
        
        tic;
        [eB,vB] = eig_sshopm(tensor(B), 'Start', xB);
        [eC,vC] = eig_sshopm(tensor(C), 'Start', xC);
        tBC(n,i) = toc;
        
        tic;
        [eA,vA] = eig_sshopm(tensor(A), 'Start', xA);
        tA(n,i) = toc;
        
        disp(string(tBC(n,i)) + " -- " + string(tA(n,i)));
        disp(sum(abs(vA - kron(vB,vC))));
    end
end

%% Make the finalized figure
set(groot, 'DefaultFigureRenderer', 'painters');
figure('Renderer', 'painters', 'Position', [0 0 900 400]);
fs = 18;

tt = tBC; tA = tA;
% Calculate mean and standard error for each row
mean_tt = mean(tt, 2);
std_error_tt = std(tt, 0, 2) / sqrt(itrs);
mean_ta = mean(tA, 2);
std_error_ta = std(tA, 0, 2) / sqrt(itrs);

% Create scatter plots
maxNP = 25;
% figure;
hold on;
errorbar(1:maxNP, mean_ta(1:maxNP), std_error_ta(1:maxNP), 'ro', 'MarkerSize', 5, 'LineWidth', 1);
errorbar(1:maxNP, mean_tt(1:maxNP), std_error_tt(1:maxNP), 'bo', 'MarkerSize', 5, 'LineWidth', 1);
hold off;

xlabel('Dimension of $$\textsf{B}$$ and $$\textsf{C}$$','Interpreter','latex');
ylabel('Time (sec.)','Interpreter','latex');
title('Z-Eigenvalue Compute Time','Interpreter','latex');
legend(["Direct", "Kronecker"]);
set(gca, 'YScale', 'log');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');
saveas(gcf, 'eigenvalueExp_08032023v1.png')