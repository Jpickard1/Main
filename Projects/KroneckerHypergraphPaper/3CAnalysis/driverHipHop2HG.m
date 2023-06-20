function driverHipHop2HG(sys, worker, epsilon, path2data, path2out)
% DRIVERHIPHOP2HG
%
%   For each hypergraph, we identify all sets of 3 vertices that are within
%   an epsilon value of one another and define this to be a multi-way
%   contact. These contacts are written as an adjacency list to a text
%   file.
%
%   EXECUTION: This is written to be executed with the corresponding bash
%   script as an array job with 15 workers. Workers 0-4 will process 'OFF',
%   5-9 'ON', and 10-14 'HIGH', then after that each worker is responsible
%   for executing the script on 20 separate files.
%
%   Usage from Joshua's laptop:
%
%       path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
%       worker = 0; epsilon = 0.01; sys = 'LT';
%       path2out = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\KroneckerHypergraphPaper\3CAnalysis\HipHop2HG\';
%
%   NOTES TO LATER ME:
%       epsilon of 0.01 appears a reasonable parameter. The number of
%       hyperedges of interest (ie not 3 adjacent vertices) appears to be
%       just over twice the number of beads
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 20, 2023

%% Preamble
disp('EXECUTE: driverHipHop2HG(sys, worker, epsilon, path2data, path2out)')
disp('         input arguments are printed out in order:');
disp(sys);
disp(worker);
disp(epsilon);
disp(path2data);
disp(path2out);
if strcmp(sys, 'GL')
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit2/'));
elseif strcmp(sys, 'DBTM')
    addpath(genpath('C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations/HyperKronFit2/'));
elseif strcmp(sys, 'LT')
    addpath(genpath('C:\Users\picka\Documents\my_rojects\DBTM\Main\Projects\AFOSR\kroneckerCalculations\HyperKronFit2\'));
else
    error(['Invalid System: the addpaths are configured for GreatLakes, the lab' ...
        'computer, or Joshuas laptop.']);
end

%% Script
% 1. Get File Path
if worker <= 4
    mid = "OFF";
elseif worker <= 9
    mid = "ON";
else
    mid = "HIGH";
end
path2data = path2data + mid + "/conf.";
suffixDNA = ".DNA";
suffixTXT = ".txt";

% 2. Set file range
fileRange = mod(worker, 5) * 40 + [1:40];

for f=1:length(fileRange)
    % Set specific file
    filePath = path2data + string(fileRange(f)) + suffixDNA;

    % Display output to user
    disp(filePath);

    % Read LAMMPS file into matlab
    T = readLAMMPSoutput(filePath);

    % Extract 
    cords = getCords(T, 5000);

    % Compute the Vietoris-Rips complex with radius parameter 0.1
    [~, F] = vrcomplex(cords, epsilon);

    % Remove all 3 way contacts between only 2 vertices
    keep = [];
    for i=1:size(F,1)
        if length(unique(F(i,:))) == 3
            keep = [keep i];
        end
    end
    F = F(keep,:);

    % Set output file path
    filePath = path2out + "adjList_" + mid + "_" + string(fileRange(f)) + "_" + string(epsilon) + suffixTXT;

    % Write output file
    writematrix(F, filePath);
end

end

