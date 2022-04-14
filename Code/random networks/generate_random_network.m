function [network] = generate_random_network(type, V, p, params)
%GENERATE_RANDOM_NETWORK: This function generates a random network. 
%   The random networks are constructed using algorithms designed to give 
%   the network specific properties that mirror real world networks. While
%   the type of network may change the network's properties, this function
%   generates networks with a specific size and number of edges, so that
%   the different kinds of networks can be compared to one another.
%
%   PARAMETERS:
%       - V: the number of vertices in the network
%       - p: the probability an edge in the network exist
%       - type: the type of random network. Current options include: 
%           - ER: Erdos-Renyi
%           - ERG: Erdos-Renyi-Gilbert
%           - SF: Scall Free
%           - SW: Small World
%           - QR: Quasi Ramanujan
%       - params: additional parameters specific to a type of random
%                 network (such as beta for SW) should be passed in params,
%                 a dictionary containing the parameter as a key-value pair

% Set the parameters for a specific network type
switch type
    case 'ER'
        % Select the number of edges from a binomial distribution
        E = binornd(nchoosek(V, 2), p);
    case 'QR'
        % Select the number of edges from a binomial distribution
        E = binornd(nchoosek(V, 2), p);
        % Calcualte the mean degree of the graph
        k = 2 * E / V;
    case 'SF'
        network = scale_free_network(V, round(0.05 * V));
    case 'SW'
        % Select the number of edges from a binomial distribution
        E = binornd(nchoosek(V, 2), p);
        % Calcualte the mean degree of the graph
        k = 2 * E / V;
    otherwise
        disp('ERROR: Enter a valid type of random network:')
end

% Construct the network
switch type
    case 'ER'
        network = erdos_renyi_network(V, E);
    case 'ERG'
        network = erdos_renyi_gilbert_network(V, 0.2);
    case 'QR'
        network = quasi_ramanujan_network(V, k, params('Beta'));
    case 'SF'
        network = scale_free_network(V, round(0.05 * V));
    case 'SW'
        network = small_world_network(V, k, params('Beta'));
end

end

