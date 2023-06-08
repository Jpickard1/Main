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

% assert(sum(sum(theta)) == 1)

n0 = size(theta, 1);
n =  n0^kronExp;

N0  = n0 * ones(1, kronExp);
MOD = N0 .^ (0:kronExp-1);
MOD = flip(MOD)';

E = zeros(numE, 2);

probDist = reshape(theta, [numel(theta), 1]);
probDist = probDist / sum(probDist);
cumDist = [0; cumsum(probDist)];

for e=1:numE

    % Drop ball through Kronecker matrix    
    falls = randsrc(1,kronExp,[1:numel(theta); reshape(theta, [1, numel(theta)])]);
    % Reindex linear indices to ij indices
    fallIJ = zeros(numel(falls), 2);
    for f=1:numel(falls)
        [fallIJ(f,1), fallIJ(f,2)] = ind2sub(size(theta), falls(f));
    end
    % disp('=======');
    % disp('linIdx');
    % disp(falls);
    % disp('ijIdk');
    % disp(fallIJ);
    % Convert to 0 indexed system
    fallIJ = fallIJ - 1;
    % Convert to kronecker exponentiated system and back to 1 indexed
    % system
    v1n = sum(fallIJ(:,1) .* MOD) + 1;
    v2n = sum(fallIJ(:,2) .* MOD) + 1;
    % disp('kronIdx');
    % disp('('+string(v1n) + ',' + string(v2n) + ')');
    % Save edge value
    E(e,:) = [v1n, v2n];
end


end

%{
    for k=1:kronExp
        r = rand;
        % linIdx = find(r>cumDist, 1, 'last');
        linIdx = find(r>cumDist, 1, 'last');
        while cumDist(linIdx) == cumDist(max([linIdx - 1, 1]));
            if linIdx == 1; break; end
            linIdx = linIdx - 1;
        end
        [i, j] = ind2sub(size(theta), linIdx);
        v1(k) = i;
        v2(k) = j;
    end
    disp([v1; v2]);
%}