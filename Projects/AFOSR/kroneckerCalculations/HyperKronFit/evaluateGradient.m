%% Evaluate Gradient
%                   
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 1, 2023
function [dldt]=evaluateGradient(A, theta)

% TOREMOVE
% T = 10;

n = size(A,1);

P = perms(1:n);
    
gradients = zeros(size(P,1), 1);
for t=1:size(P,1)
    
    p = P(t,:);
    % p = samplePermutation(A, theta)
    
    
end

end
