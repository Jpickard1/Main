%% tensor eigenvalue figure

clear; close all; clc;
load("eigencalculationRunTimes.mat");
figure('Renderer', 'painters', 'Position', [0 0 900 400]); 
fs = 18;

histogram(tA); hold on;
histogram(tB + tC);
xlabel("Time (sec.)",'interpreter','latex');
ylabel("Frequency",'interpreter','latex');
title("Maximum H-Eigenvalues of $\textsf{A}\in\mathbf{R}^{4\times 4\times 4}$ Compute Time",'interpreter','latex');
legend(["Direct", "Kronecker"],'interpreter','latex')
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');
saveas(gcf, 'fixedSizeEigenvalueCalculations_07192023_v1.png')