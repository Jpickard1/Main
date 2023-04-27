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

T = 100; n1=2;
X1 = zeros(T, n1);  X1(1,:) = rand(n1,1) - 0.5;
for t=2:T
    X1(t,:) = 0.1 * ttvk(tensor(A), X1(t-1,:)') + X1(t-1,:)';
end
figure; plot(X1)

%% Quiver of A
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
quiver3(X,Y,Z,DX,DY,DZ); hold on; title('Vector Field of A\otimesA when x_4=0');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
saveas(gcf, 'stableUnstableExample.png')

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