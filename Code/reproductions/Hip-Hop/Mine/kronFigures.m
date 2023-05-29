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
D = squareform(pdist(cords)); A = D(1:625,1:625);

[Bv, Cv] = nearestKroneckerProduct(A, [25 25 25], [25 25 25]);
Ar = kron(Bv,Cv);

figure;
subplot(2,5,[1 2 6 7]);
imagesc(A); xlabel('Beads'); ylabel('Beads'); title('Distance Map');
subplot(2,5,3); imagesc(Bv); title('Factor 1');
subplot(2,5,8); imagesc(Cv); title('Factor 2');
subplot(2,5,[4 5 9 10]); imagesc(Ar); xlabel('Beads'); ylabel('Beads'); title('Reconstructed Map');


%% 3D plot

[E, F] = vrcomplex(cords(1:625,:), 0.0075);

figure; hold on; xlabel('X'); ylabel('Y'); zlabel('Z'); title('Polymer');
plot3(cords(:,1), cords(:,2), cords(:,3), '.-');
for i = 1:size(E, 1)
    plot3(cords(E(i,:), 1), cords(E(i,:), 2), cords(E(i,:), 3), 'b', 'LineWidth', 1);
end
figure; hold on; xlabel('X'); ylabel('Y'); zlabel('Z'); title('Polymer');
plot3(cords(1:625,1), cords(1:625,2), cords(1:625,3), '.-');
for i = 1:size(E, 1)
    plot3(cords(E(i,:), 1), cords(E(i,:), 2), cords(E(i,:), 3), 'b', 'LineWidth', 1);
end


[E, F] = vrcomplex(cords(1:400,:), 0.0075);
IM = HAT.hyperedge2IM(F);
HG = Hypergraph('IM',IM)
A = HG.adjTensor;


[B, S] = tkpsvd(A, 20 * ones(1,6));
Ar = superkron(B{1,1}, B{2,1});

figure;
subplot(2,5,[1 2 6 7]); [x,y,z] = ind2sub(size(A), find(A~=0)); scatter3(x,y,z,'.');
subplot(2,5,3); [x,y,z] = ind2sub(size(B{1,1}), find(B{1,1}~=0)); scatter3(x,y,z,'.');
subplot(2,5,8); [x,y,z] = ind2sub(size(B{2,1}), find(B{2,1}~=0)); scatter3(x,y,z,'.');
subplot(2,5,[4 5 9 10]); [x,y,z] = ind2sub(size(Ar), find(Ar~=0)); scatter3(x,y,z,'.');




%% Quick experiment
%% Load data
A = cell(50,1);
for i=1:50
    disp(i)
    filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf." + string(i) + ".DNA";
    T = readLAMMPSoutput(filePath);
    cords = getCords(T);
    D = squareform(pdist(cords));
    A{i} = D(1:4900,1:4900);
end

%% Perform factorization and reconstruction
B = cell(50,1);
C = cell(50,1);
Ar = cell(50,1);
for i=1:20
    disp(i);
    [Bv, Cv] = nearestKroneckerProduct(A{i}, [70 70 70], [70 70 70]);
    B{i} = Bv;
    C{i} = Cv;
    Ar{i} = kron(Bv, Cv);
end

%% Compute distance matrices
%
%   B-B
%   B-C
%   C-C
%   A-A
%   Av-Av
%   A-Av

DAAv = zeros(20,20);
for i=1:20
    disp(i);
    for j=1:20
        DAAv(i,j) = norm(A{i} - Ar{j});
    end
end

DAA = zeros(20,20);
for i=1:20
    disp(i);
    for j=1:20
        DAAv(i,j) = norm(A{i} - A{j});
    end
end

figure; imagesc(DAAv(1:13,:))

