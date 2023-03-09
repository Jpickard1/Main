function [Xhat] = fse(A, C, systemOutput, RNG)
%FSE Full State Estimator
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 7, 2023

T  = length(systemOutput);

Kf = (lqr(A',C',eye(size(A)),eye([1 1])))';
Xhat = zeros(T,size(A,1));
Xhat(1,:) = RNG * rand(size(A,1),1);
systemOutputHat = zeros(size(systemOutput));
for t=2:T
    Xhat(t,:) = Xhat(t-1,:)' + A * Xhat(t-1,:)' + Kf * (systemOutput(t-1) - systemOutputHat(t-1));
    systemOutputHat(t) = C * Xhat(t,:)';
end

end

