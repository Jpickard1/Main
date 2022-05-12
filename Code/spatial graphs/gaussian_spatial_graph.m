function [outputArg1,outputArg2] = uniform_spatial_graph(n, k, l, d, dist, plot_on)
%UNIFORM_SPATIAL_GRAPH This function constructs a unifrom spatial graph.
%   A uniform spatial graph is a graph where vertices have a probability of
%   connecting to all other vertices within a certain length.
%
% PARAMETERS
%   - n: The number of vertices in the graph
%   - k: The average degree in the graph
%   - l: The length, given as a fraction, that the nodes can connect in
%   - d: The spatial dimension the graph is being constructed in (by
%   default d=2
%   - dist: The distance measure used to compute the bounds. By default
%   dist is euclidian (currently unimplemented)
%

%n=500;k=10;l=0.2;d=2;dist = 'euclidean';

if nargin == 3
    d = 2;
    dist = 'euclidean';
elseif nargin == 4
    dist = 'euclidean';
end

sleep_time = 0.05;

A = zeros(n,n);     % Adjacency matrix

% Create coordinates
cords = rand(n, d);
if plot_on
    figure;
    g = graph(A);
    p = plot(g);
    p.XData = cords(:,1);
    p.YData = cords(:,2);
    hold on
    title('Uniform Spatial Graph')
    pause(sleep_time);
end

pairwise_dists = squareform(pdist(cords, dist));
pairwise_dists = pairwise_dists + diag(ones(n,1) * 1000);

edges = 0;
% Add edges
vx = 1;
while edges < (k * n) / 2
    % Find which vertices are within the bounds
    vx_dists = pairwise_dists(vx,:);
    possible_vxs = find(vx_dists < l);

    % Select vx at random and add edge
    vx2 = randsample(possible_vxs, 1, false);
    A(vx, vx2) = true;
    A(vx2, vx) = true;
    edges = edges + 1;

    if plot_on
        %center = cords(vx,:);
        highlight(p, vx)
        %th = 0:pi/50:2*pi;
        %xunit = l * cos(th) + center(1);
        %yunit = l * sin(th) + center(2);
        %plot(xunit, yunit);
        %pause(sleep_time);
        highlight(p, possible_vxs)
        pause(sleep_time);
        hold off;
        g = graph(A);
        p = plot(g);
        xlim([-0.25, 1.25]);
        ylim([-0.25, 1.25]);
        p.XData = cords(:,1);
        p.YData = cords(:,2);
        title('Uniform Spatial Graph');
        hold on
        pause(sleep_time);
    end

    % Switch to the next vx in order
    vx = mod(vx, n);
    vx = vx + 1;
end
end

%{
th = 0:pi/50:2*pi;
xunit = l * cos(th) + center(1);
yunit = l * sin(th) + center(2);
plot(xunit, yunit);
xlim([-0.25, 1.25])
ylim([-0.25, 1.25])
%}

