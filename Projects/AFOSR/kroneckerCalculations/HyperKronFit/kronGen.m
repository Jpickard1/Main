function E = kronGen(theta, kronExp, numE)
% KRONGEN Generator of Kronecker graphs
%
%   theta:   initiator matrix
%   kronExp: number of times kronecker exponentiation is performed
%   numE:    number of edges to be generated
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

n0 = size(theta, 1);
n =  n0^kronExp;

N0  = n0 * ones(1, kronExp);
MOD = N0 .^ (0:kronExp-1);

E = zeros(numE, 2);

probDist = reshape(theta, [1, numel(theta)]);
probDist = probDist / sum(probDist);
cumDist = [0, cumsum(probDist)];

for e=1:numE
    % Drop ball through Kronecker matrix
    v1 = zeros(1, kronExp);
    v2 = zeros(1, kronExp);
    
    for k=1:kronExp
        r = rand;
        linIdx = find(r>cumDist, 1, 'last');
        [i, j] = ind2sub(size(theta), linIdx);
        v1(k) = i;
        v2(k) = j;
    end

    v1n = sum(v1 .* MOD);
    v2n = sum(v2 .* MOD);
    
    E(e,:) = [v1n, v2n];
end


end