%% Figure Classic Differential Equations
%
%   This figure will show a few classic differential equations, and in the
%   caption we will describe the tensor structure used to create each plot.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 23, 2023

figure('Renderer', 'painters', 'Position', [0 0 1350 400]); 
tiledlayout(1,3);
fs = 18;

% Lorenz System
sigma = 10;
rho = 28;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
m = multilinearSystem('poly',p);
[~,X] = ode45(@(t, x) m.eval(x),[0 100],[1 1 1]);
nexttile; plot3(X(:,1),X(:,2),X(:,3)); title(m.title(),'Interpreter','latex','FontSize',fs);

% Van Der Pol System
p = mvpoly('type','van der pol','mu',0.9);
m = multilinearSystem('poly',p);
[~,X] = ode45(@(t, x) m.eval(x),[0 100],[0 0.1]);
% figure; plot(X(:,1),X(:,2)); title(p.title(), 'Interpreter','latex');
nexttile; plot(X(:,1),X(:,2)); title(p.title(), 'Interpreter','latex','FontSize',fs);

% SIS System
beta = 0.15;
gamma = 0.05;
p = mvpoly('type','SIS','beta',beta,'gamma',gamma,'stoch',0);
m = multilinearSystem('poly',p);
[t,X] = ode45(@(t, x) m.eval(x),[0 100],[0.95 0.05]);
% figure; hold on;
nexttile; hold on;
plot(t,X(:,1));
plot(t,X(:,2));
title(p.title(),'Interpreter','latex','FontSize',fs);

% saveas(gcf, 'fig_classic_de_v1.png');