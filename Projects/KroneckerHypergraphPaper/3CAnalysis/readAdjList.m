function [Alist] = readAdjList(filePath, headerLines)
%READADJLIST 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 6, 2023

if nargin == 1
    headerLines = 0;
end
Alist = dlmread(filePath,'\t',headerLines);

end

