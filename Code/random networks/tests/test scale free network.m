%% Sample test
close all
V = 1000;
m0 = 5;
network = scale_free_network(V, m0);
figure; plot(sort(sum(network))); ylabel('Degree'); xlabel('Vertex'); title('Degree Abundance Plot');

%% Create networks that follow differnt power laws


V = 500;
m0 = 10;
gamma = 0.5;
network = scale_free_network(V, m0, gamma);
figure; plot(sort(sum(network))); ylabel('Degree'); xlabel('Vertex'); title('Degree Abundance Plot (gamma = ' + string(gamma) + ")");


