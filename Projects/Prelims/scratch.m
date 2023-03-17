%% Prelim Scratch
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 16, 2023

A = [1 1 0;
     1 1 1;
     0 1 1];

itrs = 4; figure;
for i=1:itrs
    subplot(1, itrs, i);
    if i == 1
        Ak = A;
    else
        Ak = kron(Ak, A);
    end
    [x, y] = find(Ak == 1);
    scatter(x, y, 'filled')
end
