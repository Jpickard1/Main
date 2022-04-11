%% DEMO RANDOM NETWORKS
%   This file generates random network by calling each function in the
%   directory so you can see a sample of each.

%% Small World Network
V = 20;
k = 4;
beta = 0.2;
[network, cords] = small_world_network(V, k, beta);
subplot(3,2,3);
gplot(network, cords', '.-');
title('Small World Network','(V=20, k=4, beta=0.2)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

% Erdos-Renyi
V = 20;
E = 40;
network = erdos_renyi_network(V, E);
subplot(3,2,1);
gplot(network, cords', '.-');
title('Erdos Renyi Network','(V=20, V=40)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

% Erdos-Renyi-Gilbert
V = 20;
p = 0.21;
network = erdos_renyi_gilbert_network(V, p);
subplot(3,2,2);
gplot(network, cords', '.-');
title('Erdos-Renyi-Gilbert Network','(V=20, p=0.21)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

% Scale Free Network
V = 20;
m0 = 5;
network = scale_free_network(V, m0);
subplot(3,2,4);
gplot(network, cords', '.-');
title('Scale Free Network','(V=20, m0=5, AB algo.)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

%% Quasi Ramanujan Networks
V = 100;
k = 10;
[network, cords] = quasi_ramanujan_network(V, k);
%subplot(3,2,5)
figure;
gplot(network, cords', '.-');
title('Quasi-Ramanujan Network','(V=20, k=4)');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
disp(sum(sum(network)));

%%
saveas(gcf, 'random network demo.png')