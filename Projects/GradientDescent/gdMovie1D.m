function [M] = gdMovie1D(f, range, theta, ftheta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 19, 2023

x = range(1):0.1:range(2);
fx = zeros(length(x), 1);
for i=1:length(x); fx(i) = f(x(i)); end

h = figure;
h.Visible = 'off';
M(length(theta)) = struct('cdata',[],'colormap',[]);
for i=1:length(theta)
    plot(x, fx,'r','LineWidth',2); hold on;
    plot(theta(1:i), 0.1+ftheta(1:i),'b','Marker','o');
    % scatter(theta(i), 0.1+f(theta(i)), 'b','filled');
    xlim(range);
    M(i) = getframe; hold off;
end

end

