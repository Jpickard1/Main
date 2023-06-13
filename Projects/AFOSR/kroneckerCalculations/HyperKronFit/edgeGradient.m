function [gradient] = edgeGradientapx(n, theta, u, v)
% Evaluates the gradient of theta with respect to a single edge
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

theta = theta / sum(sum(theta));

n0 = size(theta,1);
kronExp = log(n) / log(n0);

% edgeP = edgeProbability(n, theta, u, v);
eLL = edgeLLapx(n, theta, u, v);
noEdgeLL = log(1 - exp(eLL));

% Count the number of times an entry of theta is used
count = zeros(size(theta));
for i=1:kronExp
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
        % posGrad = (c * theta(i,j)^(c - 1));
        posGrad = c / exp(theta(i,j));

        % Set gradient if edge is not present
        % negGrad = c * (edgeP / theta(i,j)) / (1 - edgeP);
        % negGrad = -c * exp(theta(i,j)) / (1 - exp(eLL));
        % negGrad = getNoEdgeDLL(i,j,theta,u,v,kronExp);
        % negGrad = - c / (exp(theta(i,j)) * (1 - noEdgeLL));
        negGrad = getNoEdgeDLL2(theta, count, i, j, eLL);

        gradient(i,j) = posGrad - negGrad;
        % gradient(i,j) = (c / theta(i,j)) - ((k - c) / (1 - theta(i,j)));
    end
end

end