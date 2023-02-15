%% RANDOM GRAPH PROCESS

%% Set up

close all
clear all
clc

% Set vars
sleep_time = 2;
n = 5000;
layout = 'rnd';
% Set cords

switch layout
    case 'circ'
        radius = 1;
        theta=linspace(0,360,n+1); theta(end)=[];
        x=radius*cosd(theta);
        y=radius*sind(theta);
        cords = [x;y];
        clear x y radius theta
    case 'rnd'
        cords = rand(2,n);
end

%% Generate a graph with no edges

figure;
A = false(n, n);
for i=1:n
    A(i,i) = true;
end
gplot(A, cords', '.-');
set(gca,'xtick',[])
set(gca,'ytick',[])

G = graph(A);
G.XData = cords(1,:)
G.YData = cords(2,:)

%% Add in edges and sleep and update
ix = [1 1 1];
while length(ix) > 1
    pause(sleep_time/100);
    % Select one edge that doesn't yet exist
    [ix, iy] = find(~A);                    % Select missing edges in A
    new_edge = randi([1 length(ix)]);       % Randomly select one edge
    A(ix(new_edge), iy(new_edge)) = true;   % Add edge into the grap
    A(iy(new_edge), ix(new_edge)) = true;
    % Update the figure
    gplot(A, cords', '.-');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

%% Add edges, color is based on connected components
ix = [1 1 1];
while length(ix) > 1
    pause(sleep_time/100);
    % Select one edge that doesn't yet exist
    [ix, iy] = find(~A);                    % Select missing edges in A
    new_edge = randi([1 length(ix)]);       % Randomly select one edge
    % Select color
    if sum(A(ix(new_edge),:)) > 0

    end
    A(ix(new_edge), iy(new_edge)) = true;   % Add edge into the grap
    A(iy(new_edge), ix(new_edge)) = true;
    % Update the figure
    gplot(A, cords', '.-');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

