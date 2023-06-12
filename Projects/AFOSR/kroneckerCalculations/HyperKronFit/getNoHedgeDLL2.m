function DLL = getNoHedgeDLL2(theta, count, idx, eLL)
%GETNOEDGEDLL
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

DLL = -1 * count(idx{:}) * (exp(eLL) / (theta(idx{:}))) / (1 - exp(eLL));

end
