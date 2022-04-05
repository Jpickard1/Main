%% Load in a network

network = sample_network_generator(30);

%% Compute popularity matrix

popularity = edge_popularity_matrix(network);

%% Partition edges by popularity

D = 3;
[known, unknown] = popularity_partition_edges(popularity, D);

%% Plot the different partitions

% Set the coordinates to plot
G = graph(network);
p = plot(G, 'layout', 'circle');
x = p.XData;
y = p.YData;
cords = [x; y];
close;

figure;
for group=1:D
    subplot(D,1,group);
    title("Group " + string(group))
    hold on;
    gplot(known{group}, cords', '-*blue');
    gplot(unknown{group}, cords', '-red');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
end

saveas(gcf, 'popularity based edge partitioning.png');
