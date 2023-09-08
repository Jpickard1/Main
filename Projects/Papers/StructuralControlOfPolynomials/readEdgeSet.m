function [HG] = readEdgeSet(fileName)
%READEDGESET
%
%   HG = readEdgeSet('ES1ex.txt')
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

HG = struct;
HG.E = table();

fid = fopen(fileName);
tline = fgetl(fid);
numVx = -1;
while ischar(tline)
    if ~strcmp(tline(1), "#")
        if numVx ~= -1
            E = split(string(tline));
            Et= double(split(E(1), ','));
            Eh= double(split(E(2), ','));
            HG.E = [HG.E; {Et Eh}];
        else
            numVx = double(string(tline));
            HG.V = 1:numVx;
        end
    end
    tline = fgetl(fid);
end

