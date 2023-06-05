% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 24, 2023

clear; clc; close all;

A11 = [-1.2593 0.5534;
        0.5543 -0.5185];
A12 = [0.5543 -0.5185;
       -0.5185 -0.1386];
A21 = [0.5543 -0.5185;
       -0.5185 -0.1386];
A22 = [-0.5185 -0.1386;
       -0.1386 -0.7037];

A = zeros(2,2,2,2);
A(:,:,1,1) = A11;
A(:,:,1,2) = A12;
A(:,:,2,1) = A21;
A(:,:,2,2) = A22;

clearvars -except A

%% continuous trajectories
T = 100; n1=2;
s = 2;
s1 = 0.05;
s2 = 0.01;
X1 = zeros(T, n1);  X1(1,:) = s * rand(n1,1) - (s/2);
for t=2:T
    X1(t,:) = X1(t-1,:)' + s1 * ttvk(tensor(A), X1(t-1,:)');
end
X2 = zeros(T, n1^2);  X2(1,:) = kron(X1(1,:),X1(1,:)); %s * rand(n1^2,1) - (s/2);
for t=2:T
    X2(t,:) = X2(t-1,:)' + s2 * ttvk(tensor(B), X2(t-1,:)');
end

t1 = [0:(T-1)] * s1;
t2 = [0:(T-1)] * s2;

figure;
subplot(1,2,1); plot(t1, X1); ylabel('State'); xlabel('Time'); title('')
subplot(1,2,2); plot(t2, X2); ylabel('State'); xlabel('Time'); title('Kronecker System')

%% Discrete trajectories

T = 100; n1=2;
s = 1;
s1 = 0.05;
s2 = 0.01;
X1 = zeros(T, n1);  X1(1,:) = s * rand(n1,1) - (s/2);
for t=2:T
    X1(t,:) = ttvk(tensor(A), X1(t-1,:)');
end
X2 = zeros(T, n1^2);  X2(1,:) = kron(X1(1,:),X1(1,:)); %s * rand(n1^2,1) - (s/2);
for t=2:T
    X2(t,:) = ttvk(tensor(B), X2(t-1,:)');
end

t1 = [0:(T-1)] * s1;
t2 = [0:(T-1)] * s2;

figure;
subplot(1,2,1); plot(t1, X1); ylabel('State'); xlabel('Time'); title('')
subplot(1,2,2); plot(t2, X2); ylabel('State'); xlabel('Time'); title('Kronecker System')

%% Continuous
llim = -10;
ulim = 10;
spacing = 2;
[X,Y] = meshgrid(llim:spacing:ulim);
DX = zeros(size(X)); DY = zeros(size(Y));
for i=1:length(X)
    for j=1:length(X)
        v = 2 * ttvk(tensor(A), [X(i,j) Y(i,j)]');
        DX(i,j) = v(1);
        DY(i,j) = v(2);
    end
end

figure('Renderer', 'painters', 'Position', [0 0 900 400]);
subplot(1,2,1);
quiver(X,Y,DX,DY,'k','LineWidth',1.25); hold on; title('Vector Field of $\textsf{A}$','interpreter','latex');
xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex');
axis square
% contour(X,Y,Z)
% axis equal
% hold off

% Quiver of B
B = superkron(A,A);
% llim = -10;
% ulim = 10;
spacing = 4;
[X,Y,Z] = meshgrid(llim:spacing:ulim);
DX = zeros(size(X)); DY = zeros(size(Y)); DZ = zeros(size(Z));
for i=1:length(X)
    for j=1:length(X)
        for k=1:length(X)
            v = 2 * ttvk(tensor(B), [X(i,j,k) Y(i,j,k) Z(i,j,k) 0]');
            DX(i,j,k) = v(1);
            DY(i,j,k) = v(2);
            DZ(i,j,k) = v(3);
        end
    end
end
subplot(1,2,2);
quiver3(X,Y,Z,DX,DY,DZ,'k','LineWidth',1.25); hold on; title('Vector Field of $\textsf{A}\otimes\textsf{A}$ when $x_4=0$','interpreter','latex');
xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex'); zlabel('$x_3$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
% saveas(gcf, 'stableUnstableExample_05042023.png')

%% Discrete
llim = -10;
ulim = 10;
spacing = 2;
[X,Y] = meshgrid(llim:spacing:ulim);
DX = zeros(size(X)); DY = zeros(size(Y));
for i=1:length(X)
    for j=1:length(X)
        cord = [X(i,j) Y(i,j)]';
        v = 2 * ttvk(tensor(A), cord);
        DX(i,j) = cord(1) - v(1);
        DY(i,j) = cord(2) - v(2);
    end
end

figure('Renderer', 'painters', 'Position', [0 0 900 400]);
subplot(1,2,1);
quiver(X,Y,DX,DY); hold on; title('Vector Field of A');
xlabel('x_1'); ylabel('x_2');

% Quiver of B
B = superkron(A,A);
spacing = 4;
[X,Y,Z] = meshgrid(llim:spacing:ulim);
DX = zeros(size(X)); DY = zeros(size(Y)); DZ = zeros(size(Z));
for i=1:length(X)
    for j=1:length(X)
        for k=1:length(X)
            cord = [X(i,j,k) Y(i,j,k) Z(i,j,k) 0]';
            v = 2 * ttvk(tensor(B), [X(i,j,k) Y(i,j,k) Z(i,j,k) 0]');
            DX(i,j,k) = cord(1) - v(1);
            DY(i,j,k) = cord(2) - v(2);
            DZ(i,j,k) = cord(3) - v(3);
        end
    end
end
subplot(1,2,2);
quiver3(X,Y,Z,DX,DY,DZ); hold on; title('Vector Field of B=A\otimesA when x_4=0');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');

% saveas(gcf, 'stableUnstableExample.png')

%% Simpler matrix

x = [1/sqrt(2) 1/sqrt(2)];
y = [1/sqrt(2) -1/sqrt(2)];
A = x'*y;

llim = -10;
ulim = 10;
spacing = 2;
[X,Y] = meshgrid(llim:spacing:ulim);
DX = zeros(size(X)); DY = zeros(size(Y));
for i=1:length(X)
    for j=1:length(X)
        v = 2 * ttvk(tensor(A), [X(i,j) Y(i,j)]');
        DX(i,j) = v(1);
        DY(i,j) = v(2);
    end
end
figure('Renderer', 'painters', 'Position', [0 0 900 400]);
subplot(1,2,1);
quiver(X,Y,DX,DY); hold on; title('Vector Field of A');
xlabel('x_1'); ylabel('x_2');
B = superkron(A,A);
% llim = -10;
% ulim = 10;
spacing = 4;
[X,Y,Z] = meshgrid(llim:spacing:ulim);
DX = zeros(size(X)); DY = zeros(size(Y)); DZ = zeros(size(Z));
for i=1:length(X)
    for j=1:length(X)
        for k=1:length(X)
            v = 2 * ttvk(tensor(B), [X(i,j,k) Y(i,j,k) Z(i,j,k) 0]');
            DX(i,j,k) = v(1);
            DY(i,j,k) = v(2);
            DZ(i,j,k) = v(3);
        end
    end
end
subplot(1,2,2);
quiver3(X,Y,Z,DX,DY,DZ); hold on; title('Vector Field of A\otimesA when x_4=0');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');

%% superkron
B = superkron(A,A);
T = 100; n2=n1*n1;
X2 = zeros(T, n2);  X2(1,:) = rand(n2,1) - 0.5;
for t=2:T
    X2(t,:) = 0.005 * ttvk(tensor(B), X2(t-1,:)') + X2(t-1,:)';
end
figure; plot(X2);
set(ax, 'XTick', [-2 -1 0 1 2], ...
    'XTickLabel', {'-100', '-10', '0', '10', '100'});

set(gca, 'YScale', 'log')


disp(B)