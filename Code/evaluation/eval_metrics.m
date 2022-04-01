function metric = eval_metrics(known_network, unknown_network, predicted_network)
    %predicted_network = (predicted_network & ~known_network); % Predicted links only
    unknown_edges = network_to_edges_eval(unknown_network);
    predicted_edges = network_to_edges_eval((predicted_network & ~known_network));
    [~,~,~,AUC] = perfcurve(double(unknown_edges), double(predicted_edges), 1);
    metric = AUC;
end
