%% Kronecker Figures
%
%   Here I describe the structure of a few LAMMPS simulations with the
%   tensor Kronecker product
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

clear all; clc; close all;

filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);

figure;
plot3(cords(1:625,1), cords(1:625,2), cords(1:625,3), '.-');
hold on; xlabel('X'); ylabel('Y'); zlabel('Z'); title('Polymer');

n = 625;

D = squareform(pdist(cords)); A = D(1:625,1:625);   % pairwise distance matrix
D1 = repmat(A, [1, 1, n]);                          % repmat to make a 3d tensor of pairwise distances
D2 = permute(D1, [1 3 2]);                          % transpose distance matrix
D3 = permute(D1, [2 1 3]);                          % transpose distance matrix again
T = D1 + D2 + D3;                                   % get sum of pairwise distances

% Plot tensor 
