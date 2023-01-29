function [vars, multirelations] = multirelation(data, k, type)
% MULTIRELATION
%
%   This function computes the k way multirelations for all sets of k
%   variables in a data matrix.
%
% Auth: Joshua Pickard
% Date: September 14, 2022

[~, n] = size(data);

vars = nchoosek(1:n, k);
multirelations = zeros(length(vars), 1);

for i=1:length(vars)
    d = data(:, vars(i,:));
    switch type
        case 'zvi'
            multirelations(i) = zviDrezner(d);
        case 'wang-zheng'
            multirelations(i) = wangZheng(d);
        case 'taylor'
            multirelations(i) = benjaminTaylor(d);
        case 'jpic1'
            multirelations(i) = jpic1Multirelation(d);
    end
end

end
