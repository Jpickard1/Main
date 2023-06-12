function [llDiff] = hPPRtest2(p, theta, A, u, v)
% HPPRTEST2
%
% ASSUMPION:
%   1. A is a supersymmetric tensor
%
%   1. Make permutations p1 and p2
%   2. Evaluate likelihoods of rows i and j of A assuming no edges
%   3. Make updates for existing edges
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

theta = theta / sum(theta, 'all');

l1 = 0;                         % Likelihood of permutation 1
l2 = 0;                         % Likelihood of permutation 2

k = length(size(A));            % order of uniform hypergraph
n = size(A,1);                  % number of vertices in A
n0 = size(theta,1);             % size of initiator matrix
kronExp = log(n) / log(n0);     % number of repeated kron(theta, theta) until
                                % it is the size of A

theta2 = theta .^ 2;            % squares every element of theta
thetaS = zeros(n0,1);
theta2S = zeros(n0,1);
str = ",:"; for vx=1:k-2; str = str + ",:"; end
for vx=1:n0
    thetaS(vx) = eval("sum(theta(" + string(vx) + str + "),'all')");
    theta2S(vx) = eval("sum(theta2(" + string(vx) + str + "),'all')");
end

%   1. Make permutations p1 and p2
p1 = p;                 % permutation 1
p2 = p;
p2([u v]) = p2([v u]);  % permutation 2

%   2. Evaluate likelihoods of rows i and j of A assuming no edges
%       2.a Evaluate p1 first
p1u = p1(u);
p1v = p1(v);
p1uK = kronIndices(p1u, n, n0);
p1vK = kronIndices(p1v, n, n0);
llP1u = - (sum(prod(thetaS(p1uK)))) - 0.5 * (sum(prod(theta2S(p1uK))));
llP1v = - (sum(prod(thetaS(p1vK)))) - 0.5 * (sum(prod(theta2S(p1vK))));
l1 = l1 + llP1u + llP1v;

%       2.a Evaluate p2 first
p2u = p2(u);
p2v = p2(v);
p2uK = kronIndices(p2u, n, n0);
p2vK = kronIndices(p2v, n, n0);
llP2u = - (sum(prod(theta(p2uK)))) - 0.5 * (sum(prod(theta2S(p2uK))));
llP2v = - (sum(prod(theta(p2vK)))) - 0.5 * (sum(prod(theta2S(p2vK))));
l2 = l2 + llP2u + llP2v;

%   3. Make updates for existing edges
E = getHyperedges(A);
E = sort(E,2);
E = unique(E, 'rows');
hContainsUV = find(sum((E == v) + (E == u), 2) > 0);

for i=1:length(hContainsUV)
    idx = E(hContainsUV(i), :);
    p1idx = num2cell(p1(idx));
    p2idx = num2cell(p2(idx));

    eLLp1 = hedgeLLapx(n, theta, p1idx);
    eLLp2 = hedgeLLapx(n, theta, p2idx);

    ll1update = -log(1-exp(eLLp1)) + eLLp1;
    ll2update = -log(1-exp(eLLp2)) + eLLp2;

    l1 = l1 + ll1update;
    l2 = l2 + ll2update;
end

llDiff = l2 - l1;

end
