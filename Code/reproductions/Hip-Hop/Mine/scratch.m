% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

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


