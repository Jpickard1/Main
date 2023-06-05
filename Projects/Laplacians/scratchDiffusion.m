%% Scratch Diffusion
%
%   This file is meant to simulate and visualize diffusion using graph and
%   hypergraph laplacians.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 4, 2023

clear all;
close all;
clc;
pon = true;

%% Diffusion along a graph chain with n components
n = 20;
A = diag(ones(n-1,1),1); A = A + A'; A(n,1) = 1; A(1,n) = 1;
D = diag(sum(A));
L = D - A;

T = 10;
% X0 = 10 * rand(n,1);
X0 = zeros(n,1); X0(1) = 100;
[t, X] = ode45(@(t,y)diffusionODE(t,y,L), 1:T, X0,L);
if pon; figure; plot(X); title('Agreement/Diffusion Protocol'); xlabel('Time'); ylabel('Node State'); end
if pon; figure; plot(graph(A)); title('Chain Graph'); end

%% Make a movie
M(T) = struct('cdata',[],'colormap',[]);
h = figure;
% h.Visible = 'off';
for t=1:T
    % scatter(1:n, X(:,t));
    % ylim([0, 10.05]);
    p = plot(graph(A));
    for v=1:n
       highlight(p,v,'NodeColor', [1 1 1-(X(t,v) / sum(X(t,:)))]);
    end
    drawnow;
    M(t) = getframe;
end
% h.Visible = 'on';
figure; title('Movie'); movie(M);

%% ode45
T = 1000;
X = zeros(n, T);
X0 = rand(n,1)';
[t, X] = ode45(@(t,y)diffusionODE(t,y,L), 1:T, X0,L);

%% Hypergraph Diffusion
HG = hyperchain(n,3);
L = HG.laplacianTensor;
X0 = 2 * rand(n,1)';
[t, X] = ode45(@(t,y)diffusionHG(t,y,L), 1:T, X0,L);
if pon; figure; plot(X); title('Hypergraph: Agreement/Diffusion Protocol'); xlabel('Time'); ylabel('Node State'); end
A = diag(ones(n-1,1),1); A = A + A'; D = diag(sum(A)); L = D - A;
[t, X] = ode45(@(t,y)diffusionODE(t,y,L), 1:T, X0,L);
if pon; figure; plot(X); title('Graph: Agreement/Diffusion Protocol'); xlabel('Time'); ylabel('Node State'); end


function dy = diffusionODE(t, y, L)
    dy = - L * y;
end

function dy = diffusionHG(t, y, L)
    dy = -1 * ttvk2(tensor(L), y, 1);
end
