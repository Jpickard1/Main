%% Experiment Parameters

graph_type = ["ER", "ERG", "SW", "SF"];
edge_removal_method = ["R", "SB", "CE", "RC"];
V = 100;

%% Experiment

for gt=graph_type
    for er=edge_removal_method
        network = generate_random_network(gt, 100);
        figure;
        subplot(3,1,1);
        scree(double(network));
        title(string(gt) + " True");
        [known, unknown] = edge_removal(network, er, 0.5);
        subplot(3,1,2);
        scree(double(known));
        title(string(er) + " Removal: Known");
        subplot(3,1,3);
        scree(double(unknown));
        title(string(er) + " Removal: Unknown");
    end
end


