%% Figure for Kronecker Hypergraph Paper
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 20, 2023

clear; close all; clc
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.2.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

% Compute the Vietoris-Rips complex with radius parameter 0.1
[E, F] = vrcomplex(cords, 0.01);
FF = F;
keep = [];
for i=1:size(F,1)
    if length(unique(F(i,:))) == 3
        keep = [keep i];
    end
end
F = F(keep,:);
keep = find(range(F,2) > 7);
F = F(keep,:);

figure;
subplot(1,3,1);
plot3(cords(:,1), cords(:,2), cords(:,3), '.-','LineWidth',1,'MarkerSize',3);
title('Pax6 Polymer Structure'); axis off;
subplot(1,3,2);
plot3(cords(:,1), cords(:,2), cords(:,3), '.-','LineWidth',1,'MarkerSize',3); hold on;
for i = 1:size(F, 1)
    fill3(cords(F(i,:), 1), cords(F(i,:), 2), cords(F(i,:), 3), 'r', 'EdgeColor','none'); hold on;
end
title('Chromatin Multi-way Contacts'); axis off;
subplot(1,3,3);
% plot3(cords(:,1), cords(:,2), cords(:,3), '.-','LineWidth',1,'MarkerSize',3); hold on;
for i = 1:size(F, 1)
    fill3(cords(F(i,:), 1), cords(F(i,:), 2), cords(F(i,:), 3), 'r', 'EdgeColor','none'); hold on;
end
title('Observed Multi-way Contacts'); axis off;
saveas(gcf, 'polymer_and_triangles.png');

%%
IM = HAT.hyperedges2IM(F);
HG = Hypergraph('IM', IM);
A = HG.cliqueGraph;

cc = conncomp(graph(A));
vxs = find(cc == mode(cc));
keep = find(sum(ismember(F,vxs),2) > 0);

F = F(keep,:);
IM = HAT.hyperedges2IM(F);
HG = Hypergraph('IM', IM);

%%
figure; HG.plot()

figure; HG.plot();

subplot(1,3,3);
for i = 1:size(F, 1)
    fill3(cords(F(i,:), 1), cords(F(i,:), 2), cords(F(i,:), 3), 'r', 'EdgeColor','none'); hold on;
end
title('Highlighted Multi-way Contacts'); axis off;

% xlabel('X'); ylabel('Y'); zlabel('Z'); title('Polymer');

