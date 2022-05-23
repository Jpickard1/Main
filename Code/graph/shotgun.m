function [outputArg1,outputArg2] = shotgun(A)
%SHOTGUN Call this function first for any data matrix. It produces 3 plots
%   you should look at.
%
% Auth: Joshua Pickard jpic@umich.edu
% Date: May 23, 2022

figure;
subplot(2,2,1); scree(A);
subplot(2,2,2); plot(sort(sum(A), 'descend'));
%subplot(2,2,3:4); plot(graph(A));
%subplot(2,2,4);

[v,e] = eigs(A, 2, 'smallestabs');
figure; plot(v(:,1))

end

gplot(A, xy)
