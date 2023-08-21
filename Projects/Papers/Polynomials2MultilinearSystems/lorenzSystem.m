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
figure; plot3(X(:,1),X(:,2),X(:,3)); title(p.title(),'Interpreter','latex');

%% Classic Lorenz System (multilinearSystem)
sigma = 10;
rho = 28;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
m = multilinearSystem('poly',p);
[~,X] = ode45(@(t, x) m.eval(x),[0 100],[1 1 1]);
figure; plot3(X(:,1),X(:,2),X(:,3)); title(m.title(),'Interpreter','latex');

%% What to show in case studies?
%   - how do eigenvalues change at bifurcation points
%       fix beta, vary sigma and rho, show bifurcation point

deltaS = 0.5;
deltaR = 0.1;

beta = 8/3;
Svals = 5:deltaS:15;
Rvals = 22:deltaR:25;

E = cell(length(Svals), length(Rvals));

for si = 1:length(Svals)
    sigma = Svals(si);
    for ri = 1:length(Rvals)
        rho = Rvals(ri);
        disp(string(sigma) + " -- " + string(rho));
        p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
        % [~,X] = ode45(@(t, x) p.eval(x),[0 100],[1 1 1]);
        % figure; plot3(X(:,1),X(:,2),X(:,3)); title(m.title(),'Interpreter','latex');
        m = multilinearSystem('poly',p);
        A = double(m.A);
        evals = zeig(A);
        E{si, ri} = evals;
    end
end

Emax = zeros(size(E));
for si = 1:length(Svals)
    for ri = 1:length(Rvals)
        Emax(si,ri) = max(E{si, ri});
    end
end

figure; imagesc(Emax)

%% Eigenvalues of A as a sigma changes

E = {};
sigma = 10;
% rho = 28;
beta = 8/3;

for rho = 23.25:0.001:24
    disp(rho);
    p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
    % [~,X] = ode45(@(t, x) p.eval(x),[0 100],[1 1 1]);
    % figure; plot3(X(:,1),X(:,2),X(:,3)); title(m.title(),'Interpreter','latex');
    m = multilinearSystem('poly',p);
    A = double(m.A);
    evals = zeig(A);
    E{length(E)+1} = evals;
end

EE = [];
for i=1:length(E)
    EE = [EE; max(E{i})];
end


figure; plot(23.25:0.001:24, EE)


hA = hosvd(tensor(A),1e-4)


hosvd(A,1e-4)




