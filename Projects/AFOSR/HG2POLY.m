function [systemStr] = HG2POLY(HG, vars)
%HG2POLY Creates polynomial system of hypergraph
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 1, 2023

if nargin == 1
    vars = repmat(["x"], size(IM, 1), 1);
    for i=1:length(vars); vars(i) = vars(i) + "_" + string(i); end;
end

IM = HG.IM;
systemStr = "";
for vx=1:size(IM,1)
    vxDir = "\dot{" + vars(vx) + "} & = ";
    incidentHyperedges = find(IM(vx,:) ~= 0);
    for e=1:length(incidentHyperedges)
        hyperedge = find(IM(:,incidentHyperedges(e)) ~= 0);
        hyperedge = hyperedge(hyperedge ~= vx);
        if e ~= 1
            vxDir = vxDir + " + ";
        end
        for i=1:length(hyperedge)
            vxDir = vxDir + vars(hyperedge(i));
        end
    end
    systemStr = systemStr + vxDir;
    if vx ~= size(IM, 1); systemStr = systemStr + "\\" + newline; end;
end
end