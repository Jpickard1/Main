%% Wound Healing Scratch
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 28, 2023


%% Satureday 04/01/2023
clear; close all; clc;
v = VideoReader('2023-03-24-triangle-stitched.avi');
f = read(v,1);
img = f;

v.NumFrames

% F = kron(f, ones(e,e));
e=4;
D = rgb2gray(f); D = double(D); D = kron(D, ones(e,e));
figure; imagesc(D); hold on;

BW = edge(D)

Dt = (D > 0.3 * mean(mean(f)));

figure; imagesc(f); hold on;
[x, y] = find(Dt == 0);
sc = scatter(y, x, '.');

figure; imagesc(f)




%% Tuesday 
% vidObj = VideoReader('youtubeWoundClip.mp4', 'CurrentTime',1.2);

imdata = imread('Screenshot 2023-03-28 164047.jpeg');
D = rgb2gray(imdata);
Dt = (D > 0.8 * mean(mean(D)));
figure; imagesc(imdata); hold on;
[x, y] = find(Dt == 0);
sc = scatter(y, x, '.');

%% Find islands from Dt
clear;
imdata = imread('Screenshot 2023-03-28 164047.jpeg');
% imdata = imread('Screenshot 2023-03-28 180412.jpeg');
D = rgb2gray(imdata);
Dt = (D > 0.8 * mean(mean(D)));
figure; imagesc(D)

[AG,NSV,G] = myIslands(Dt);
AG(AG == 47) = 0;
figure; imagesc(imdata); hold on;
[x, y] = find(AG ~= 0);
sc = scatter(y, x, '.');

cords = [];
uqs = unique(AG);
for i=1:length(uqs)
    if uqs(i) == 0; continue; end;
    [x, y] = find(AG == uqs(i));
    cords(i, 1) = mean(x);
    cords(i, 2) = mean(y);
end

scatter(cords(:, 2), cords(:, 1), '.');

%%

newIsland = 1;
islands = zeros(size(Dt));
for i=1:size(islands, 1)
    for j=1:size(islands, 2)
        if Dt(i, j) == 1
            if i ~= 1 && j ~= 1 && i < size(islands, 1)
                islandNum = max([islands(i-1, j-1), islands(i-1, j), islands(i, j-1), islands(i+1, j-1)]);
            elseif i~=1
                islandNum = islands(i-1, j);
            elseif j~=1
                islandNum = islands(i, j-1);
            else % If they are both 1
                islandNum = 1;
            end
            if islandNum == 0
                islandNum = newIsland;
                newIsland = newIsland + 1;
            end
            islands(i, j) = islandNum;
        end
    end
end
unq = unique(islands);
for i=1:length(unq)
    [x, y] = find(islands == unq(i));
    if numel(x) < 4
        islands(x, y) = 0;
    end
end

figure; imagesc(imdata); hold on;
[x, y] = find(islands,  1);
sc = scatter(y, x, '.');


%%
alpha(sc,.2)

sc.MarkerFaceAlpha = .2;




image([xmin xmax], [ymin ymax],data)

figure; imagesc(imdata)

D = rgb2gray(imdata);

figure; imagesc(D)

mean(mean(D))

Dt = (D > 0.8 * mean(mean(D)));
figure; imagesc(Dt);

figure; imagesc(imdata); hold on;
[x, y] = find(Dt == 0);
scatter(x, y, '.');
