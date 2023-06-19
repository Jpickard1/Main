function [llDiff] = permutationProbabilityRatioTest(p, idxs, E, theta)
%PERMUTATIONPROBABILITYRATIOTEST 


l1 = 0;                         % Likelihood of permutation 1
l2 = 0;                         % Likelihood of permutation 2

k = size(E,2);                   % k is order of the hypergraph
n = max(max(E));                  % number of vertices in A
n0 = size(theta,1);             % size of initiator matrix
% kronExp = log(n) / log(n0);     % number of repeated kron(theta, theta) until
                                % it is the size of A
theta2 = theta .^ 2;            % squares every element of theta

if k==2                         % Track if our hypergraph is directed or 
    directed = true;            % undirected

    thetaR = sum(theta, 2);         % row sum of theta
    thetaC = sum(theta, 1);         % column sum of theta
    theta2R = sum(theta2, 2);       % row sum of theta squared
    theta2C = sum(theta2, 1);       % colum sum of theta squared

    thetaSC = [thetaR'; thetaC];
    theta2SC = [theta2R'; theta2C];
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
p1 = p;
p2 = p;
p2([idxs(1) idxs(2)]) = p2([idxs(2) idxs(1)]);

%   2. Evaluate likelihoods of rows idxs(1) and idxs(2) of A assuming no edges
p1u = p1(idxs(1));
p1v = p1(idxs(2));
p2u = p2(idxs(1));
p2v = p2(idxs(2));
p1uK = kronIndices(p1u, n, n0);
p1vK = kronIndices(p1v, n, n0);
p2uK = kronIndices(p2u, n, n0);
p2vK = kronIndices(p2v, n, n0);
for i=1:size(thetaSC,1)
    tS  = thetaSC(i,:);
    t2S = theta2SC(i,:);
    
    llP1u = - (sum(prod(tS(p1uK)))) - 0.5 * (sum(prod(t2S(p1uK))));
    llP1v = - (sum(prod(tS(p1vK)))) - 0.5 * (sum(prod(t2S(p1vK))));
    l1 = l1 + llP1u + llP1v;
    
    llP2u = - (sum(prod(tS(p2uK)))) - 0.5 * (sum(prod(t2S(p2uK))));
    llP2v = - (sum(prod(tS(p2vK)))) - 0.5 * (sum(prod(t2S(p2vK))));
    l2 = l2 + llP2u + llP2v;
end
% If the graph is directed, then duplicate edges were counted, so we undo
% this effect here
if directed
    l1 = l1 + log(1 - exp(hedgeLL(n, theta, p1(idxs)))) + ...
              log(1 - exp(hedgeLL(n, theta, p1(idxs)))) + ...
              log(1 - exp(hedgeLL(n, theta, p1(idxs)))) + ...
              log(1 - exp(hedgeLL(n, theta, p1(idxs))));
    l2 = l2 + log(1 - exp(hedgeLL(n, theta, p2(idxs)))) + ...
              log(1 - exp(hedgeLL(n, theta, p2(idxs)))) + ...
              log(1 - exp(hedgeLL(n, theta, p2(idxs)))) + ...
              log(1 - exp(hedgeLL(n, theta, p2(idxs))));
end

%   3. Make updates for existing edges
eIdxs = find(sum(ismember(E,idxs), 2));
for i=1:size(eIdxs)
    hedge = E(eIdxs(i), :);

    eLLp1 = hedgeLL(n, theta, p1(hedge));
    eLLp2 = hedgeLL(n, theta, p2(hedge));

    ll1update = - log(1-exp(eLLp1)) + eLLp1;
    ll2update = - log(1-exp(eLLp2)) + eLLp2;

    l1 = l1 + ll1update;
    l2 = l2 + ll2update;
end

%   4. Compute difference in log likelihoods
llDiff = l2 - l1;

end

