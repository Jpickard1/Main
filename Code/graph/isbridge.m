function [bridge] = isbridge(network, i, j)
%ISBRIDGE This function returns if the edge between nodes i and j is a
%   bridge in the network. This code assumes the network is undirected and
%   fully connected.

bridge = true;
initial_connected_components = conncomp(graph(network));
network(i,j) = false;
network(j,i) = false;
if initial_connected_components == conncomp(graph(network))
    bridge = false;
end

end
