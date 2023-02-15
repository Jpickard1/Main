function HGo = observedHG(HG)
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

V = find(sum(full(HG.IM),2) > 0);
E = find(sum(full(HG.IM),1) > 1);
HGo = Hypergraph('IM', HG.IM(V,E));

end