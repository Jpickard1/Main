function [network] = generate_random_network(type, V, arg1)
%GENERATE_RANDOM_NETWORK: This function generates a random network of the 
%   specefied type. The number of nodes in the network is set by V, and the
%   following are valid types of networks:
%   - ER: Erdos-Renyi
%   - ERG: Erdos-Renyi-Gilbert
%   - LW: Large World
%   - SF: Scall Free
%   - SW: Small World
switch type
    case 'ER'
        network = erdos_renyi_network(V, round(10 * V));
    case 'ERG'
        network = erdos_renyi_gilbert_network(V, 0.2);
    case 'LW'
        disp('LW needs to be implemented');
    case 'SF'
        network = scale_free_network(V, round(0.05 * V));
    case 'SW'
        network = small_world_network(V, 5);
    otherwise
        disp('ERROR: Enter a valid type of random network:')
end

end

