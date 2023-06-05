%% HyperKronFit
%              
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [theta]=HyperKronFit(A)

lr = 0.01;

n = size(A,1);
k = length(size(A));

theta = rand(2 * ones(1, k));

maxItrs = 10;
for itr=1:maxItrs
    dldt = evaluateGradient(theta, A);
    theta = theta + lr * dldt;
end

end
