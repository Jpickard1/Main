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
        network = erdos_renyi_network(V, arg1);
    case 'ERG'
        network = erdos_renyi_gilbert_network(V, arg2);
    case 'LW'
        disp('LW needs to be implemented');
    case 'SF'
        network = scale_fee_network(V);
    case 'SW'
        disp('SW needs to be implemented');
    otherwise
        disp('ERROR: Enter a valid type of random network:')
end

end

