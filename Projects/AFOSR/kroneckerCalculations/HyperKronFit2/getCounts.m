function [count] = getCounts(n, theta, idxs)
%GETCOUTNS Counts the contribution of each element in theta to a Kronecker 
% expanded index.

n0 = size(theta, 1);
kronExp = ceil(log(n) / log(n0));

% Count the number of times an entry of theta is used
count = zeros(size(theta));
c = cell(kronExp, 1);
kIdx = zeros(kronExp, length(idxs));
for i=1:length(idxs)
    kIdxs(:,i) = kronIndices(idxs(i), n, n0);
end
for i=1:kronExp
    c = num2cell(kIdxs(i,:));
    count(c{:}) = count(c{:}) + 1;
end

% for i=1:kronExp
%     i1 = mod(floor((idxs(1) - 1)/ n0^(i-1)), n0) + 1;
%     i2 = mod(floor((idxs(2) - 1)/ n0^(i-1)), n0) + 1;
%     count(i1, i2) = count(i1, i2) + 1;
% end

end

