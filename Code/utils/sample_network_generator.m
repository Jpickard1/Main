function network = sample_network_generator(number_of_nodes)
    %network = randi([0 1], number_of_nodes, number_of_nodes);
    adj = rand(number_of_nodes);
    adj = (adj > 0.5);
    network = tril(adj) | tril(adj)';
    for i=1:number_of_nodes
        network(i,i) = false;
    end
end
