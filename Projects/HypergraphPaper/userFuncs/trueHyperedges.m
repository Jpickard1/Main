function [U] = trueHyperedges(HGtrue, Uidx)
%TRUEHYPEREDGES Extract true hyperedges for hyperedge prediction
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: January 29, 2023

U = cell(length(Uidx), 1);
for i=1:length(U)
    e = find(HGtrue.IM(:, Uidx(i)) ~= 0);
    if length(e) < 2
        continue
    end
    U{i} = e;
end


end

