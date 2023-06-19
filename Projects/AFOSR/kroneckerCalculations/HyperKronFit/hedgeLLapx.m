function [LL] = hedgeLLapx(n, theta, idxs)
%HEDGELLAPX Approximates the log likelihood of a hyperedge
%
% NOTE: This function works for an arbitraty k-uniform hypergraph, so I
% should always use this function for graphs as well.
%
% PARAMETERS
%   n is number of vertices in a hypergraph
%   theta is the kronecker initiator parameters
%   idxs are the vertices of the hyperedge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 12, 2023

if ismatrix(theta)
    LL = edgeLL(n, theta, idxs(1), idxs(2));
    return;
end

n0 = size(theta,1);
kronExp = log(n) / log(n0);

LL = 0;
idx = zeros(kronExp, numel(idxs));
for i=1:kronExp
    idx(i,:) = mod(floor((idxs - 1) / n0^(i - 1)), n0) + 1;
end
if ndims(theta) == 2
    linIdx = sub2ind(size(theta), idx(:,1),idx(:,2));
elseif ndims(theta) == 3
    linIdx = sub2ind(size(theta), idx(:,1),idx(:,2));
else
    error('Joshua, you wrote a bad function here');
end
LL = LL - sum(theta(linIdx)) - 0.5 * sum(theta(linIdx).^2);

% for i=1:kronExp
%     idx = mod(floor((idxs - 1) / n0^(i - 1)), n0) + 1;
%     idx = num2cell(idx);
%     LL = LL - theta(idx{:}) - 0.5 * theta(idx{:})^2;
% end


end