%% Gradient Descent Scratch File
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 19, 2023

f = @(x) x^2;
df = @(x) 2 * x;
x = -2:0.4:2;
fx = zeros(length(x),1); for i=1:length(x); fx(i) = f(x(i)); end

learningRate = 0.9;

itrs = 20;
theta = zeros(itrs,1); ftheta = zeros(size(theta));
theta(1) = 2; ftheta(1) = f(theta(1));
for i=2:itrs
    g = df(theta(i - 1));
    theta(i) = theta(i-1) - learningRate * g;
    ftheta(i) = f(theta(i));
end

M = gdMovie1D(f, [-2, 2], theta, ftheta);
figure; movie(M)

figure; hold on;
plot(x, fx,'r','LineWidth',2,'Marker','o','');
figure; scatter(theta, ftheta);

%% 3D gradient descent
clear; close all; clc;
% x = linspace(-2*pi,2*pi);
% y = linspace(-2*pi,2*pi);
x = 0:0.01:4*pi;
y = x;
[X,Y] = meshgrid(x,y);
% Z = sin(X) + cos(Y);
% contour(X,Y,Z);

[X, Y, Z] = peaks;
learningRate = 0.7;
itrs = 20;
theta = zeros(itrs,2); theta(1,:) = randsample(x, 2);
thetax = zeros(itrs,2);
thetay = zeros(itrs,2);
thetaz = zeros(itrs,2);
for i=2:itrs
    g = [cos(theta(i-1, 1)) -sin(theta(i-1,2))];
    theta(i,:) = theta(i-1,:) - learningRate * g;
    [~, xi] = find(X(1,:) < theta(i,1), 1, 'last');
    [~, yi] = find(Y(:,1) < theta(i,2), 1, 'last');
    thetaz(i) = Z(xi, yi);
    thetax(i) = xi; thetay(i) = yi;
end

figure;
contour3(X,Y,Z,20); hold on;
plot(theta(:,1), theta(:,2), 'b')

figure; surf(X, Y, Z); hold on;
scatter3(thetay, thetax, thetaz)

% scatter(theta(:,1), theta(:,2));