function [llDiff] = PPRtest3(p, theta, A, u, v, E)
% PPRTEST3
%
%   This function performes a permutation probability ratio test for either
%   graphs or hypergraphs. It extends the previous PPRtest2.m function that
%   worked only for graphs.
%
%   This function evaluates equations 5.7 and 5.8 in the thesis to
%   determine if the likelihood of the adjacency matrix being generated by
%   A will increase or decrease when the labeling of As nodes is evaluated
%   with permutations p and p with the terms u and v swapped.
%
%   1. Make permutations p1 and p2
%   2. Evaluate likelihoods of rows i and j of A assuming no edges
%   3. Make updates for existing edges
%   4. Compute difference in log likelihoods
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 15, 2023

if nargin == 4
    E = getEdgesFromAdj(A);
end

l1 = 0;                         % Likelihood of permutation 1
l2 = 0;                         % Likelihood of permutation 2

k = ndims(A);                   % k is order of the hypergraph
n = size(A,1);                  % number of vertices in A
n0 = size(theta,1);             % size of initiator matrix
kronExp = log(n) / log(n0);     % number of repeated kron(theta, theta) until
                                % it is the size of A
theta2 = theta .^ 2;            % squares every element of theta

if k==2                         % Track if our hypergraph is directed or 
    directed = true;            % undirected

    thetaR = sum(theta, 2);         % row sum of theta
    thetaC = sum(theta, 1);         % column sum of theta
    theta2R = sum(theta2, 2);       % row sum of theta squared
    theta2C = sum(theta2, 1);       % colum sum of theta squared

    thetaSC = {thetaR thetaC};
    theta2SC = {theta2R theta2C};
else
    directed = false;

    thetaS = zeros(n0,1);
    theta2S = zeros(n0,1);
    str = ",:"; for vx=1:k-2; str = str + ",:"; end
    for vx=1:n0
        thetaS(vx) = eval("sum(theta(" + string(vx) + str + "),'all')");
        theta2S(vx) = eval("sum(theta2(" + string(vx) + str + "),'all')");
    end

    thetaSC = {thetaS};
    theta2SC = {theta2S};
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
for i=1:size(thetaSC,2)
    tS = thetaSC{i};
    t2S = theta2SC{i};
    llP1u = - (sum(prod(tS(p1uK)))) - 0.5 * (sum(prod(t2S(p1uK))));
    llP1v = - (sum(prod(tS(p1vK)))) - 0.5 * (sum(prod(t2S(p1vK))));
    l1 = l1 + llP1u + llP1v;
end

%       2.a Evaluate p2 first
p2u = p2(u);
p2v = p2(v);
p2uK = kronIndices(p2u, n, n0);
p2vK = kronIndices(p2v, n, n0);
for i=1:size(thetaSC,2)
    tS = thetaSC{i};
    t2S = theta2SC{i};
    llP2u = - (sum(prod(tS(p2uK)))) - 0.5 * (sum(prod(t2S(p2uK))));
    llP2v = - (sum(prod(tS(p2vK)))) - 0.5 * (sum(prod(t2S(p2vK))));
    l2 = l2 + llP2u + llP2v;
end

% If the graph is directed, then duplicate edges were counted, so we undo
% this effect here
if directed
    l1 = l1 + log(1 - exp(edgeLLapx(n, theta, p1(u), p1(u)))) + ...
              log(1 - exp(edgeLLapx(n, theta, p1(u), p1(v)))) + ...
              log(1 - exp(edgeLLapx(n, theta, p1(v), p1(u)))) + ...
              log(1 - exp(edgeLLapx(n, theta, p1(v), p1(v))));
    l2 = l2 + log(1 - exp(edgeLLapx(n, theta, p2(u), p2(u)))) + ...
              log(1 - exp(edgeLLapx(n, theta, p2(u), p2(v)))) + ...
              log(1 - exp(edgeLLapx(n, theta, p2(v), p2(u)))) + ...
              log(1 - exp(edgeLLapx(n, theta, p2(v), p2(v))));
end

%   3. Make updates for existing edges
% E = getEdgesFromAdj(A);
idxs = find(sum((E == u) + (E == v), 2) > 0);
for i=1:size(idxs)
    edge = E(idxs(i), :);

    eLLp1 = hedgeLLapx(n, theta, p1(edge));
    eLLp2 = hedgeLLapx(n, theta, p2(edge));

    ll1update = -log(1-exp(eLLp1)) + eLLp1;
    ll2update = -log(1-exp(eLLp2)) + eLLp2;

    l1 = l1 + ll1update;
    l2 = l2 + ll2update;
end

%   4. Compute difference in log likelihoods
llDiff = l2 - l1;


end