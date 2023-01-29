function [multirelationTensor] = multirelationTensor(data, k, type)
% MULTIRELATION TENSOR 
%
%   This function constructs a k-way multirelation tensor from a data
%   matrix.
%
% Auth: Joshua Pickard
% Date: September 14, 2022

[~, n] = size(data);

% Compute multirelation
[idxs, mrlns] = multirelation(data, k, type);

% Create tensor
T = tenzeros(n*ones(1,k));
multirelationTensor = symtensor(T);

% Fill tensor
for i=1:length(idxs)
    multirelationTensor(idxs(i,:)) = mrlns(i);
end

end

