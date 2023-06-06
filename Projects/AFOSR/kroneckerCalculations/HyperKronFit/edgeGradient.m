function [gradient] = edgeGradient(n, theta, u, v)
% Evaluates the gradient of theta with respect to a single edge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

n0 = 2;

k = log(n) / log(n0);

p = 1;

% Count the number of times an entry of theta is used
count = zeros(size(theta));
for i=1:k
    i1 = mod(floor((u - 1)/ n0^(i-1)), n0) + 1;
    i2 = mod(floor((v - 1)/ n0^(i-1)), n0) + 1;
    count(i1, i2) = count(i1, i2) + 1;
end

gradient = zeros(size(theta));
% Calculate gradient of log likelihood function
for i=1:size(theta,1)
    for j=1:size(theta,2)
        % Count the number of times (i,j) was used
        c = count(i,j);

        % Set gradient if edge is present
        posGrad = (c * theta(i,j)^(c - 1));
        
        % Set gradient if edge is not present
        posGrad = (c * theta(i,j)^(c - 1));

        gradient(i,j) = gradient(i,j) + posGrad - negGrad;
        % gradient(i,j) = (c / theta(i,j)) - ((k - c) / (1 - theta(i,j)));
    end
end

end