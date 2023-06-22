function [E] = loadDAWN(sys, k)
%LOADDAWN Constructs adjacency list from DAWN dataset
%
% sys: which computer is the program running on
% k  : order of the hypergraph
%
%   sys = "DBTM"; k = 3;
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 21, 2023

if strcmp(sys, 'DBTM')
    path2data = "C:\Joshua\MissingData\Projects\KroneckerHypergraphPaper\BigFit\DAWN\";
    path2out = "C:\Joshua\MissingData\Projects\KroneckerHypergraphPaper\BigFit\adjLists\";
elseif strcmp(sys, 'GL')
    path2data = "/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/BigFit/DAWN/";
    path2out = "/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/BigFit/adjLists/";
else
    error('Invalid system');
end

nverts    = readtable(path2data + "DAWN-nverts.txt");
simplices = readtable(path2data + "DAWN-simplices.txt");

nverts = nverts{:,:};
simplices = simplices{:,:};

E = zeros(length(simplices), k);

vxCounter = 1;
numHyperedges = 1;
for i=1:size(nverts,1)
    simplexSize = nverts(i);
    if simplexSize < k
        vxCounter = vxCounter + simplexSize;
        continue;
    end
    vertices = simplices(vxCounter:vxCounter+simplexSize-1);
    simplexKclique = nchoosek(vertices,k);
    E(numHyperedges:numHyperedges + size(simplexKclique,1)-1,:) = simplexKclique;
    vxCounter = vxCounter + simplexSize;
    numHyperedges = numHyperedges + size(simplexKclique,1);
end

outfile = path2out + "DAWN_" + string(k) + ".txt";
writematrix(E, outfile);

end