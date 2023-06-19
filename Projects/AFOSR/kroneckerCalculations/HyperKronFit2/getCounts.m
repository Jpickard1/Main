function [count] = getCounts(n, theta, idxs)
%GETCOUTNS Counts the contribution of each element in theta to a Kronecker 
% expanded index.

n0 = size(theta, 1);
kronExp = log(n) / log(n0);

% Count the number of times an entry of theta is used
count = zeros(size(theta));
for i=1:kronExp
    i1 = mod(floor((idxs(1) - 1)/ n0^(i-1)), n0) + 1;
    i2 = mod(floor((idxs(2) - 1)/ n0^(i-1)), n0) + 1;
    count(i1, i2) = count(i1, i2) + 1;
end

end

