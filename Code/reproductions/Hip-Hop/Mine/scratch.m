% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

%% April 12, 2023 vrcomplex and LAMMPS output
clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));
neighborDistances = diag(D,1); max(neighborDistances)

(sum(sum(D <= 0.0075)) - 5000) / 2

% Compute the Vietoris-Rips complex with radius parameter 0.1
[E, F] = vrcomplex(cords, 0.0075);

figure; hold on; xlabel('X'); ylabel('Y'); zlabel('Z'); title('Polymer');
plot3(cords(:,1), cords(:,2), cords(:,3), '.-');
for i = 1:size(E, 1)
    plot3(cords(E(i,:), 1), cords(E(i,:), 2), cords(E(i,:), 3), 'b', 'LineWidth', 1);
end
% subplot(1,3,1); plot3(cords(:,1), cords(:,2), cords(:,3), '.-'); title('Polymer'); xlabel('X'); ylabel('Y'); zlabel('Z');
% subplot(1,3,2); imagesc(D); title('Distance Matrix'); xlabel('Loci'); ylabel('Loci');
% subplot(1,3,3); imagesc(HC); title('Pseudo Hi-C Matrix'); xlabel('Loci'); ylabel('Loci');

%% April 12, 2023 vrcomplex verification
clear; close all; clc

% Generate a 2D point cloud with 100 points
points = rand(100,2);

% Compute the Vietoris-Rips complex with radius parameter 0.1
[E, F] = vrcomplex(points, 0.1);

% Plot the point cloud and edges
figure;
scatter(points(:,1), points(:,2), 'k', 'filled');
hold on;
for i = 1:size(E, 1)
    plot(points(E(i,:), 1), points(E(i,:), 2), 'b', 'LineWidth', 1);
end
axis equal;
title('Vietoris-Rips complex with r=0.1');

% Plot the triangles
% figure;
% scatter(points(:,1), points(:,2), 'k', 'filled');
hold on;
for i = 1:size(F, 1)
    patch(points(F(i,:), 1), points(F(i,:), 2), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
axis equal;
title('Triangles in Vietoris-Rips complex with r=0.1');



%% April 12, 2023 Pseudo Pore-C

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));
P = (D < mean(mean(D)));
%%

epsilon = mean(mean(D));
n = size(D,1);
S = [];
for i=1:n
    for j=i+1:n
        % If loci i and j are too far apart, then no simplicial complices
        % can contain them
        if P(i,j) == 0; continue; end
        neighbors = unique(intersect(find(P(i,:) == 1), find(P(j,:) == 1)));
        T = zeros(length(neighbors), 3); T(:,1) = i; T(:,2) = j;
        T(:,3) = neighbors;
        S = [S; T];
        if mod(j,1000) == 0; disp(j); end
    end
    disp(i)
end

%%
figure; 
subplot(1,3,1); plot3(cords(:,1), cords(:,2), cords(:,3), '.-'); title('Polymer'); xlabel('X'); ylabel('Y'); zlabel('Z');
subplot(1,3,2); imagesc(D); title('Distance Matrix'); xlabel('Loci'); ylabel('Loci');
subplot(1,3,3); imagesc(HC); title('Pseudo Hi-C Matrix'); xlabel('Loci'); ylabel('Loci');

%% April 11, 2023 Pseudo Pore-C

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));
HC = (D < mean(mean(D)));

figure; 
subplot(1,3,1); plot3(cords(:,1), cords(:,2), cords(:,3), '.-'); title('Polymer'); xlabel('X'); ylabel('Y'); zlabel('Z');
subplot(1,3,2); imagesc(D); title('Distance Matrix'); xlabel('Loci'); ylabel('Loci');
subplot(1,3,3); imagesc(HC); title('Pseudo Hi-C Matrix'); xlabel('Loci'); ylabel('Loci');


%% Find all 3 way simplices
n = size(D,1);
for i=1:n
    for j=i:n
        % If loci i and j are too far apart, then no simplicial complices
        % can contain them
        if HC(i,j) == 0; continue; end
        neighbors = unique(union(HC(i,:), HC(j,:)));
        for ki=1:length(neighbors)
            k = neighbors(ki);
            vi = cords(i,:);
            vj = cords(j,:);
            vk = cords(k,:);
            mean = mean([vi; vj; vk])

        end
    end
end

n = 1000;
G = cell(n, n, n);
X = [min(cords(:,1)) max(cords(:,1))];
Y = [min(cords(:,2)) max(cords(:,2))];
Z = [min(cords(:,3)) max(cords(:,3))];
X = X(1):(X(2)-X(1))/n:X(2);
Y = Y(1):(Y(2)-Y(1))/n:Y(2);
Z = Z(1):(Z(2)-Z(1))/n:Z(2);



%% April 11, 2023 PLOTTING

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

figure; 
subplot(1,2,1); plot3(cords(:,1), cords(:,2), cords(:,3), '.-'); title('Polymer'); xlabel('X'); ylabel('Y'); zlabel('Z');
subplot(1,2,2); imagesc(D); title('Distance Matrix'); xlabel('Loci'); ylabel('Loci');

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


