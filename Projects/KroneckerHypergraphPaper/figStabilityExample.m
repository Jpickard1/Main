%% Stability Example

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
rng(1)
T = 10000; n1=2;
s = 1;
s1 = 0.075;
s2 = 0.075;
AA = superkron(A,A); % + 0.01 * rand(4,4,4);

t1 = [0:(T-1)] * s1;
t2 = [0:(T-1)] * s2;

% System
X1 = zeros(T, n1);  X1(1,:) = s * rand(n1,1) - (s/2);
for t=2:T
    X1(t,:) = X1(t-1,:)' + s1 * ttvk(tensor(A), X1(t-1,:)');
end
% Kronecker system
X2 = zeros(T, n1^2);  X2(1,:) = kron(X1(1,:),X1(1,:)); %s * rand(n1^2,1) - (s/2);
for t=2:T
    X2(t,:) = X2(t-1,:)' + s2 * ttvk(tensor(AA), X2(t-1,:)');
end

figure;
tiledlayout(1,2);
nexttile;
plot(t1, X1); ylabel('State'); xlabel('Time'); title('Base System');
nexttile;
plot(t2, X2); ylabel('State'); xlabel('Time'); title('Kronecker System'); 
set(gca, 'YScale', 'log');

% % Kronecker system + noise
% X3 = zeros(T, n1^2);  X2(1,:) = kron(X1(1,:),X1(1,:)); %s * rand(n1^2,1) - (s/2);
% for t=2:T
%     X3(t,:) = X3(t-1,:)' + s2 * ttvk(tensor(AA), X3(t-1,:)');
% end
% % picking of adequate values for the labels
% TickMask = linspace(1,numel(y),N);
% YTickLabels = y(TickMask);
% % scale labels and plotdata, remove NaN ->inconsistency, do you really want that?
% YTick = scale( YTickLabels );
% Y = scale(y);
% YTick(isnan(YTick)) = 0;
% Y(isnan(Y)) = 0;
% set(gca,'YTick',YTick,'YTickLabels',YTickLabels)
%% 
set(groot, 'DefaultFigureRenderer', 'painters');
figure('Renderer', 'painters', 'Position', [0 0 900 900]);
tiledlayout(2,2);
fs = 18;

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

nexttile;
quiver(X,Y,DX,DY); hold on; title('Vector Field of $\dot{\mathbf{x}}=\textsf{B}\mathbf{x}^3$','interpreter','latex');
xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex');
axis square
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');

% contour(X,Y,Z)
% axis equal
% hold off

% Quiver of B
B = superkron(A,A);
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
nexttile;
quiver3(X,Y,Z,DX,DY,DZ); hold on; title('Vector Field of $\dot{\mathbf{x}}=\textsf{A}\mathbf{x}^3$ when $x_4=0$','interpreter','latex');
xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex'); zlabel('$x_3$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');

nexttile;
plot(t1, X1); ylabel('State'); xlabel('Time'); title('Example Trajecytory of $\dot{\mathbf{x}}=\textsf{B}\mathbf{x}^3$','interpreter','latex');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');

nexttile;
plot(t2, X2); ylabel('State'); xlabel('Time'); title('Example Trajectory of $\dot{\mathbf{x}}=\textsf{A}\mathbf{x}^3$','interpreter','latex'); 
set(gca, 'YScale', 'log');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');
saveas(gcf, 'stableUnstableExample_07192023.png')
