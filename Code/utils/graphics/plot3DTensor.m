function [outputArg1,outputArg2] = plot3DTensor(T)
%PLOT3DTENSOR Summary of this function goes here
%   Detailed explanation goes here
%
% Auth: Joshua Pickard
% Date: May 28, 2023


n = size(T,1);
cords = nchoosek(1:n,3);
c = num2cell(cords',2);
D = [cords T(sub2ind(size(T),c{:}))'];

% figure;
scatter3(D(:,1),D(:,2),D(:,3),40,D(:,4),'filled')    % draw the scatter plot

end

