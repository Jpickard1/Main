function [cords] = getCords(T, maxBead)
%GETCORDS Get coordinates of beads of interest from LAMMPS output struct
%
%   This function accepts a struct generated from readLAMMPSoutput and
%   returns the 3D coordinates of the beads of interests. Importan
%
% EXAMPLE:
%   T     = readLAMMPSoutput(filePath);
%   cords = getCords(T);
%
% Auth: Joshua B. Pickard
%       jpic@umich.edu
% Date: April 9, 2023

if nargin == 1; maxBead = 5000; end;

idxs = (T.id <= maxBead);
cords = [T.xs T.ys T.zs];   % Coordinates of all beads
cords = cords(idxs,:);      % Extract only beads along the string
[~, I] = sort(T.id(idxs));  % Sort beads so that the cords is ordered by
cords = cords(I,:);         % the bead number

end

