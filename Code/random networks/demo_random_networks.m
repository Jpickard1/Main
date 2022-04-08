%% DEMO RANDOM NETWORKS
%   This file generates random network by calling each function in the
%   directory so you can see a sample of each.

%% Small World Network
V = 20;
k = 4;
beta = 0.2;
[network, cords] = small_world_network(V, k, beta);
subplot(2,2,3);
gplot(network, cords', '.-');
title('Small World Network','(V=20, k=4, beta=0.2)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

% Erdos-Renyi
V = 20;
E = 40;
network = erdos_renyi_network(V, E);
subplot(2,2,1);
gplot(network, cords', '.-');
title('Erdos Renyi Network','(V=20, V=40)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

% Erdos-Renyi-Gilbert
V = 20;
p = 0.21;
network = erdos_renyi_gilbert_network(V, p);
subplot(2,2,2);
gplot(network, cords', '.-');
title('Erdos-Renyi-Gilbert Network','(V=20, p=0.21)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

% Scale Free Network
V = 20;
m0 = 5;
network = scale_free_network(V, m0);
subplot(2,2,4);
gplot(network, cords', '.-');
title('Scale Free Network','(V=20, m0=5, AB algo.)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

saveas(gcf, 'random network demo.png')