function [HGp] = observeV(HG, r, mp)
%OBSERVER
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

IM = full(HG.IM);
[n, m] = size(IM);

% Initialize vertex and hyperedge lists (line 2 in tex)
Eprime = [];
Vprime = [];

% Bias and normalize popularity (line 3 in tex)
p = sum(IM,2);
p = p .^ r;
p(isinf(p)) = 0;
p = p / sum(p);

hyperedgeCardinality = sum(IM,1);

% (line 4 in tex)
while length(Eprime) < mp
    % (line 5 in tex)
    i = randsample(1:n,1,true,p);
    % (line 6 in tex)
    Vprime = [Vprime i];
    % (line 7 in tex)
    Eprime = find(sum(IM(Vprime,:), 1) >= 0.9 * hyperedgeCardinality);
    % (line 8 in tex)
    p(i) = 0;
    p = p / sum(p);
    %disp(sum(p));
    %disp(string(length(Vprime)) + ": " + string(length(Eprime)));
end

% Remove excess edges
%{
while length(Eprime) >= mp
    E = find(IM(Vprime(end),:) ~= 0);
    removeE = randsample(E,1);
    Eprime = Eprime(Eprime ~= removeE);
end
%}

% Extract hypergraph (lines 10 in tex)
IMp = zeros(n,m);
IMp(Vprime,:) = IM(Vprime,:);
% IMp(:,Eprime) = IM(:,Eprime);

% Remove hyperedges incident to only one vertex (line 7 in tex)
% singleVxEdges = find(sum(IMp) == 1);
% IMp(:,singleVxEdges) = 0;

HGp = Hypergraph('IM',IMp);

if size(observedHG(HGp).IM,2) > mp
    disp('Catch')
end

end

%% Testing
%{
HG = HAT.uniformErdosRenyi(5,4,3)
IM = full(HG.IM)
r = 1
mp = 3
%}