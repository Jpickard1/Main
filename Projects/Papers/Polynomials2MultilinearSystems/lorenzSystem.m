%% Lorenz System
%
%   This file constructs the lorenz system as a multilinear system.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 7, 2023

%% Classic Lorenz System (mvpoly)
sigma = 10;
rho = 2;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
[~,X] = ode45(@(t, x) p.eval(x),[0 1000],[1 1 1]);
figure; plot3(X(:,1),X(:,2),X(:,3)); title(p.title(),'Interpreter','latex');

%% Classic Lorenz System (multilinearSystem)
sigma = 10;
rho = 28;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
m = multilinearSystem('poly',p);
[~,X] = ode45(@(t, x) m.eval(x),[0 100],[1 1 1]);
figure; plot3(X(:,1),X(:,2),X(:,3)); title(m.title(),'Interpreter','latex');

%% August 21, 2023

sigma = 10;
rho = 28;
beta = 8/3;
D = [];
for rho=20:0.1:28
    p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
    m = multilinearSystem('poly',p);
    A = double(m.A);
    A = tensor(A);
    C = cp_als(A, 1);
    D = [D; C.lambda * [C.U{1}; C.U{2}; C.U{3}]'];
end

Y = tsne(D)

species = [ones(41,1); zeros(40,1)];
figure;
gscatter(Y(:,1),Y(:,2), species,eye(3))
title('2-D Embedding')

%% Lorenze System Hopf and Pitchfork Bifurcations

S = [];
D = [];
for rho=0.7:0.1:1.5
    disp("rho: " + string(rho))
    for sigma=8:0.2:12
        disp("sigma: " + string(sigma))
        for beta=8/3:0.5:15
            % disp(beta)
            p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
            m = multilinearSystem('poly',p);
            A = double(m.A);
            A = tensor(A);
            C = cp_als(A, 1,'printitn',0);
            D = [D; C.lambda * [C.U{1}; C.U{2}; C.U{3}]'];

            if rho < 1
                S = [S; 0];
            elseif sigma > beta + 1
                S = [S; 1];
            else
                S = [S; 2];
            end
        end
    end
end

Y = tsne(D);
species = S;
% species = [ones(41,1); zeros(40,1)];
% species = [3 * ones(1721,1); S];
figure;
gscatter(Y(:,1),Y(:,2), species,eye(3))
title('2-D Embedding')

C = corr(D');
C(isnan(C)) = 0;

[V, ~] = eigs(C, 2);

figure;
gscatter(V(:,1),V(:,2), species,eye(3))
title('Spectral Embedding')

figure; imagesc(C)

%% Chaotic bifurcations
%   this bifurcates near rho=23.5

sigma = 10;
beta = 8/3;
S = [];
D = [];

for rho = 22.0:0.005:25
    disp(rho);
    p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
    % [~,X] = ode45(@(t, x) p.eval(x),[0 100],[1 1 1]);
    % figure; plot3(X(:,1),X(:,2),X(:,3)); title(m.title(),'Interpreter','latex');
    m = multilinearSystem('poly',p);
    A = double(m.A);
    A = tensor(A);
    C = cp_als(A, 1,'printitn',0,'tol',1e-8);
    D = [D; [C.U{1}; C.U{2}; C.U{3}]'];
    if rho < 23.5
        S = [S; 0];
    else
        S = [S; 1];
    end
end

% spectral embeddings
C = corr(D');
C(isnan(C)) = 0;

[V, ~] = eigs(C, 2);

figure;
plot(V(:,1),V(:,2)); %, species,eye(2))
title('Spectral Embedding')

figure; imagesc(C)

% tSNE visualization
Y = tsne(D);
species = S;
% species = [ones(41,1); zeros(40,1)];
% species = [3 * ones(1721,1); S];
figure;
gscatter(Y(:,1),Y(:,2), species,eye(3))
title('tSNE Embedding')

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




