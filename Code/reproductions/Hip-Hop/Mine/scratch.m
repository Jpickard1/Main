% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

%% Pseudo Pore-C
%   I would like to say k points are in a contact if all k points fit
%   within a sphere of radius epsilon.
%
%   As a starting point (k=3), I have the pairwise distances between all
%   points. 3 points are within a sphere, or circle in this case, when...
%
%   How are topological bar codes computed? That could be used here
%
%   I could do a grid search over the 3D space of each bead. For each bead,
%   I could select all grid coordinates within the beads bounds at a given
%   epsilon, and then I could mark that grid as being contained in that
%   bead. Then, I could find intersections by looking at the individual
%   boxes of my grid.
%
%       Grid: min(cords) to max(cords)
%       n   : number of boxes in each dimension
%       eps : epsilon value

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

n = 100;
G = cell(n, n, n);

X = [min(cords(:,1)) max(cords(:,1))];
Y = [min(cords(:,2)) max(cords(:,2))];
Z = [min(cords(:,3)) max(cords(:,3))];
X = X(1):(X(2)-X(1))/n:X(2);
Y = Y(1):(Y(2)-Y(1))/n:Y(2);
Z = Z(1):(Z(2)-Z(1))/n:Z(2);

for bi=1:size(cords,1)
    % Select bead
    b = cords(bi,:);
    % Find box containing the bead
    xi = find(X < b(1), 1, 'last'); % xi is X coordinate of grid box
    zi = find(Y < b(2), 1, 'last');
    yi = find(Z < b(3), 1, 'last');

    % Mark bead b as containing xi zi and yi in G
    G{xi, zi, yi} = [G{xi, zi, yi} bi];

    g = [X(xi) Y(yi) Z(zi)];

    norm(b-g)
    
end


%% Pseudo Hi-C
%   ISSUE: the read in LAMMPS dump file does not contain the beads numbered
%   along the string.

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);

figure; plot3(cords(:,1), cords(:,2), cords(:,3))

D = squareform(pdist(cords));
figure; imagesc(D); title('Simulated Exact Distances');

A = D(1:4900,1:4900);
[B, C] = nearestKroneckerProduct(D, [2 2], [2500 2500]);
figure; subplot(2,1,1); imagesc(B); subplot(2,1,2); imagesc(C);

thresh = mean(mean(D));
phc = (D > thresh);
figure; imagesc(phc); title('Pseudo Hi-C')

%% Test readLAMMPSoutput

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = [T.xs T.ys T.zs];
idxs = (T.id <= 5000);

figure; hold on;
for id=1:max(T.type(:))
    ids = find(T.type(idxs,:) == id);
    plot3(cords(ids,1),cords(ids,2),cords(ids,3), '.')
end


figure; scatter3(cords(idxs,1),cords(idxs,2),cords(idxs,3), '.')


