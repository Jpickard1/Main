function cords = circleCords(A)
%CIRCLECORDS This function assigns coordinates to plotting the graph given
%   by A. When A is a scalar, the coordinates are assigned randomly. When A
%   is a matrix, the coordinates are assigned so that the vertices are
%   plotted in degree order.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2022

% Select number of vertices
n = length(A);
if length(A) == 1
    n = A;
end

% Set coordinates in circle
radius = 1;
theta=linspace(0,360,n+1); theta(end)=[];
x=radius*cosd(theta);
y=radius*sind(theta);
cords = [x;y]';

figure; plot(graph(A), 'XData', cords(:,1), 'YData', cords(:,2))
% Assign coordinates in order of degree
if length(A) ~= 1
    degree = sum(A);
    [degree, sortIdx] = sort(degree);
    % sort B using the sorting index
    cords = cords(sortIdx, :);
    
    figure; plot(graph(A), 'XData', cords(:,1), 'YData', cords(:,2))
end

end

