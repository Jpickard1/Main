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

theta = theta / sum(theta, 'all');

n0 = size(theta, 1);
n =  n0^kronExp;
k = ndims(theta);

N0  = n0 * ones(1, kronExp);
MOD = N0 .^ (0:kronExp-1);
MOD = flip(MOD)';

E = zeros(numE, k);

probDist = reshape(theta, [numel(theta), 1]);
probDist = probDist / sum(probDist);
cumDist = [0; cumsum(probDist)];

for e=1:numE

    % Drop ball through Kronecker matrix    
    falls = randsrc(1,kronExp,[1:numel(theta); reshape(theta, [1, numel(theta)])]);
    % Reindex linear indices to ij indices
    idxs = cell(1, k);
    [idxs{:}] = ind2sub(n0 * ones(1, k), falls');
    fallIJ = cell2mat(idxs);
    % Convert to 0 indexed system
    fallIJ = fallIJ - 1;
    % Convert to kronecker exponentiated system and back to 1 indexed
    % system
    hedge = zeros(k,1);
    for vx=1:k
        hedge(vx) = sum(fallIJ(:,vx) .* MOD) + 1;
    end
    % Save edge value
    E(e,:) = hedge;
end

end

% k = length(size(theta));
% n0 = size(theta, 1);
% n =  n0^kronExp;
% 
% N0  = n0 * ones(1, kronExp);
% MOD = N0 .^ (0:kronExp-1);
% 
% E = zeros(numE, k);
% 
% probDist = reshape(theta, [1, numel(theta)]);
% probDist = probDist / sum(probDist);
% cumDist = [0, cumsum(probDist)];
% 
% for e=1:numE
% 
%     % Drop ball through Kronecker matrix
%     vk = zeros(k, kronExp);
% 
%     for kItr=1:kronExp
%         r = rand;
%         indices = cell(k, 1);
%         linIdx = find(r>cumDist, 1, 'last');
%         [indices{1:k}] = ind2sub(size(theta), linIdx);
%         vk(:,kItr) = cell2mat(indices);
%     end
% 
%     hyperedge = zeros(1,k);
%     for vx=1:k
%         hyperedge(vx) = sum(vk(k,:) .* MOD);
%     end
% 
%     E(e,:) = hyperedge;
% end