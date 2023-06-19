function [idx0] = kronIndices(idx, n, n0)
%KRONINDICES Maps indices of a matrix
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 11, 2023

kronExp = ceil(log(n) / log(n0));
idx0 = zeros(1, kronExp);
for i=1:kronExp
    idx0(i) = mod(floor((idx - 1)/ n0^(i-1)), n0) + 1;
end
idx0 = flip(idx0);

end

%% A few tests
%{

idx = 3; n = 4; n0 = 2;
kronIndices(idx,n,n0)

%}