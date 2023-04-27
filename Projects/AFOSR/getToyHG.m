function [HG] = getToyHG(n, k, t)
%GETTOYHG This function returns a toy hypergraph
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 2023
    switch t
        case "hyperring"
            HG = hyperring(n, k);
        case "hyperchain"
            HG = hyperchain(n, k);
        case "hyperstar"
            HG = hyperstar(n, k);
        case "complete"
            HG = completeHG(n, k);
    end
end

