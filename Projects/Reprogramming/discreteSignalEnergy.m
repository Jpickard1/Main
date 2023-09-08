function [E] = discreteSignalEnergy(Y)
%DISCRETESIGNALENERGY
%
% INPUT: Y is a p x t matrix where each column contains p measurements and
% samples are taken at t times.
%
% Reference: Learning perturbation-inducible cell states from observability
% analysis of transcriptome dynamics, Equation (12)
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

E = 0;
for i=1:size(Y,2)
    E = E + Y(:,i)' * Y(:,i);
end

end

