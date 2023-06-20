function [T] = readLAMMPSoutput(filePath)
%READLAMMPSOUTPUT Summary of this function goes here
%   Detailed explanation goes here
%
% EXAMPLE:
%   filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\ON\conf.1.DNA";
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

%% Function Body
% Open file
fid = fopen(filePath);

% Number of of timesteps
n = fgetl(fid);
tSteps = str2num(fgetl(fid));
n = fgetl(fid);
nElts = str2num(fgetl(fid));
n = fgetl(fid);
xlims = str2num(fgetl(fid));
ylims = str2num(fgetl(fid));
zlims = str2num(fgetl(fid));
n = fgetl(fid);
clear n;

beads = zeros(nElts, 8); %cell2table(cell(0,9),'VariableNames',{'ATOMS', 'id', 'type', 'xs', 'ys', 'zs', 'ix', 'iy', 'iz'});
for elt=1:nElts
    beadInfo = fgetl(fid);
    c = split(beadInfo); c = c(1:8)';
    c = cellfun(@str2num,c);
    beads(elt,:) = c;
end

% Create table
T = array2table(beads);
T.Properties.VariableNames = {'id', 'type', 'xs', 'ys', 'zs', 'ix', 'iy', 'iz'};

% Close file
fclose(fid);
end

