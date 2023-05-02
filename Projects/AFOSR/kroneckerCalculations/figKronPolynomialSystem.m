%% Kronecker Polynomial System
%
%   This figure is intended to describe the Kronecker product of polynomial
%   systems
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 1, 2023

clear all; clc; close all;

IM1 = [1 1 1 0;
      0 1 1 1];
IM1 = IM1';
H1 = Hypergraph('IM', IM1);
A1 = H1.adjTensor;
Vx = ["x_1","x_2","x_3","x_4"];

IM2 = [1 1 1];
IM2 = IM2';
H2 = Hypergraph('IM', IM2);
A2 = H2.adjTensor;
Vy = ["y_1","y_2","y_3"];

Ak = superkron(A1, A2);
IMk = HAT.A32IM(Ak);
Hk = Hypergraph('IM', IMk);
Vk = kronString(Vx, Vy);

figure;
subplot(1,3,1); HAT.plotIncidenceMatrix(H1, 'sort', false); yticks(1:4); yticklabels(Vx);
subplot(1,3,2); HAT.plotIncidenceMatrix(H2, 'sort', false); yticks(1:3); yticklabels(Vy); 
subplot(1,3,3); HAT.plotIncidenceMatrix(Hk, 'sort', false); yticks(1:12); yticklabels(Vk);


Ak = superkron(A2, A1);
IMk = HAT.A32IM(Ak);
Hk = Hypergraph('IM', IMk);
Vk = kronString(Vy, Vx);

figure;
subplot(1,3,1); HAT.plotIncidenceMatrix(H1, 'sort', false); yticks(1:4); yticklabels(Vx);
subplot(1,3,2); HAT.plotIncidenceMatrix(H2, 'sort', false); yticks(1:3); yticklabels(Vy); 
subplot(1,3,3); HAT.plotIncidenceMatrix(Hk, 'sort', false); yticks(1:12); yticklabels(Vk);

S1 = HG2POLY(H1, Vx);
S2 = HG2POLY(H2, Vy);
Sk = HG2POLY(Hk, Vk);

%%

clear all; clc; close all;

IM1 = [1 1 1];
IM1 = IM1';
H1 = Hypergraph('IM', IM1);
A1 = H1.adjTensor;
Vx = ["x_1","x_2","x_3"];

IM2 = [1 1 1];
IM2 = IM2';
H2 = Hypergraph('IM', IM2);
A2 = H2.adjTensor;
Vy = ["y_1","y_2","y_3"];

Ak = superkron(A1, A2);
IMk = HAT.A32IM(Ak);
Hk = Hypergraph('IM', IMk);
Vk = kronString(Vx, Vy);

figure;
subplot(1,3,1); HAT.plotIncidenceMatrix(H1, 'sort', false); yticks(1:3); yticklabels(Vx);
subplot(1,3,2); HAT.plotIncidenceMatrix(H2, 'sort', false); yticks(1:3); yticklabels(Vy); 
subplot(1,3,3); HAT.plotIncidenceMatrix(Hk, 'sort', false); yticks(1:9); yticklabels(Vk);

S1 = HG2POLY(H1, Vx);
S2 = HG2POLY(H2, Vy);
Sk = HG2POLY(Hk, Vk);

%% For latex figure
theta = 0:(2*pi)/12:2*pi;
x = 1.5 * cos(theta); y = 1.5 * sin(theta);
figure; scatter(x, y);


