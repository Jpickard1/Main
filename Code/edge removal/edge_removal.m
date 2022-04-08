function [known_network, unknown_network] = edge_removal(network, method, percent_removed)
    if ~islogical(network)
        network = logical(network);
        disp('WARNING: Edge removal directory is set to work with undirected, unweighted graphs.')
    end
    switch method
        case 'R'
            [known_network, unknown_network] = random_removal(network, percent_removed);
        case 'SB'
            [known_network, unknown_network] = snowball_removal(network, percent_removed);
        case 'CE'
            disp('Cold Ends method not yet impimented')
        case 'RC'
            [known_network, unknown_network] = right_censoring_removal(network, percent_removed);
    end
end