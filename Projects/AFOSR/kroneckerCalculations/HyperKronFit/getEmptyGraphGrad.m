function emptyGradient = getEmptyGraphGrad(n, theta)
%GETEMPTYGRAPHGRAD
%
%   SNAP: kronecker.cpp 1147-1155 GetApxEmptyGraphDLL
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 6, 2023

% theta = theta / sum(sum(theta));

% emptyGradient = zeros(size(theta));

n0 = size(theta,1);

kronExp = log(n) / log(n0);

s  = sum(sum(theta));
sq = sum(sum(theta .^ 2));

emptyGradient = -kronExp * (s ^ (kronExp - 1)) - kronExp * (sq ^ (kronExp - 1)) * theta;
% loop on parameters of theta
% for i=1:n0
%     for j=1:n0
%         %Calculate gradient of theta(i,j)
%         emptyGradient(i,j) = -kronExp * (s ^ (kronExp - 1)) - kronExp * (sq ^ (kronExp - 1)) * theta(i,j);
%     end
% end

end