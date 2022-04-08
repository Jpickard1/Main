%% Experiment Parameters

graph_type = ["ER", "ERG", "SW", "SF"];
edge_removal = ["R", "SB", "CE", "RC"];
V = 100;

%% Experiment

for gt=graph_type
    for er=edge_removal
        network = generate_random_network(gt, 100);
        subplot(3,1,1);
        scree(double(network));
        [known, unknown] = edge_removal(network, er, 0.2);
        scree(known);
        subplot(3,1,2);
        scree(unknown);
        subplot(3,1,3);
    end
end
