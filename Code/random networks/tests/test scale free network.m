%% Sample test
close all
V = 1000;
m0 = 5;
network = scale_free_network(V, m0);
figure; plot(sort(sum(network))); ylabel('Degree'); xlabel('Vertex'); title('Degree Abundance Plot');

%% Create networks that follow differnt power laws


V = 5000;
m0 = 10;
gamma = 0;
network = scale_free_network(V, m0, gamma);
figure; plot(sort(sum(network))); ylabel('Degree'); xlabel('Vertex'); title('Degree Abundance Plot (gamma = ' + string(gamma) + ")");

%%

V = 500;
p = 0.0064;
network = scale_free_network(V, m0, gamma);
figure; plot(sort(sum(network))); ylabel('Degree'); xlabel('Vertex'); title('Degree Abundance Plot (gamma = ' + string(gamma) + ")");

%%

figure;
leg = [];
for gamma = [0 0.01 0.1 0.25 0.5 1];
    t = 1:1000;
    p = t.^-gamma;
    hold on
    plot(p);
    leg = [leg string(gamma)];
end
legend(leg);


%%
V = 500;
m0 = 10;
figure; hold on;
leg = [];
for gamma = [0]
    network = scale_free_network(V, m0, gamma);
    D = sort(sum(network));
    plot(D);
    leg = [leg string(gamma)];
end
legend(leg);

%% 

V = 500;
p = 0.3;
network = erdos_renyi_gilbert_network(V, p);
D = sort(sum(network));
figure; plot(D);


