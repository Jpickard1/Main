%% Lorenz System
%
%   This file constructs the lorenz system as a multilinear system.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 7, 2023

%% Classic Lorenz System (mvpoly)
sigma = 10;
rho = 28;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
[~,X] = ode45(@(t, x) p.eval(x),[0 100],[1 1 1]);
figure; plot3(X(:,1),X(:,2),X(:,3)); title(p.title(),'Interpreter','latex')

%% Classic Lorenz System (multilinearSystem)
sigma = 10;
rho = 28;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
m = multilinearSystem('poly',p);
[~,a] = ode45(@(t, x) m.eval(x),[0 100],[1 1 1]);
figure; plot3(a(:,1),a(:,2),a(:,3))