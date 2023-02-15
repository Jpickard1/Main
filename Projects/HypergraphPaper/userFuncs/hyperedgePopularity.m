function [pop2] = hyperedgePopularity(IM)
%HYPEREDGEPOPULARITY Summary of this function goes here
%   Detailed explanation goes here
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: January 27, 2023

[v,e] = size(IM);
vD = sum(IM,2);
pop2 = zeros(1,e);
for i=1:e
    vxc = find(IM(:,i) ~= 0);
    pop2(i) = prod(vD(vxc)) / length(vxc);
end

% When no vertices are incident set to 0. By default MATLAB set these to 1 
% with the above loop.
eC = sum(IM, 1);
pop2(find(eC == 0)) = 0;

end

