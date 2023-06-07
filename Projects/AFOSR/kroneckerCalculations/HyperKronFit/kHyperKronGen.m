function [E] = kHyperKronGen(theta, kronExp, numE)
%HYPERKRONGEN Summary of this function goes here
%
%   theta:   initiator k mode tensor
%   kronExp: number of times kronecker exponentiation is performed
%   numE:    number of edges to be generated
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

k = length(size(theta));
n0 = size(theta, 1);
n =  n0^kronExp;

N0  = n0 * ones(1, kronExp);
MOD = N0 .^ (0:kronExp-1);

E = zeros(numE, k);

probDist = reshape(theta, [1, numel(theta)]);
probDist = probDist / sum(probDist);
cumDist = [0, cumsum(probDist)];

for e=1:numE

    % Drop ball through Kronecker matrix
    vk = zeros(k, kronExp);
    
    for kItr=1:kronExp
        r = rand;
        indices = cell(k, 1);
        linIdx = find(r>cumDist, 1, 'last');
        [indices{1:k}] = ind2sub(size(theta), linIdx);
        vk(:,kItr) = cell2mat(indices);
    end

    hyperedge = zeros(1,k);
    for vx=1:k
        hyperedge(vx) = sum(vk(k,:) .* MOD);
    end
    
    E(e,:) = hyperedge;
end

end

