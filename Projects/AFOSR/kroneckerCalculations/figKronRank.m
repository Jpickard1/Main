%% Kronecker Rank
%
%   I try to generate plots that describe the set ofpossible rank one
%   matrices.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 29, 2023

%% Kronecker decomposed matrices

pts = 100000;
data = zeros(pts, 4);
for i=1:pts
    x = rand(2,1);
    y = rand(2,1);
    data(i,:) = kron(x,y);
end

figure;
labels = ["x_1","x_2","x_3","x_4"];
[h,ax] = plotmatrix(data);                        % create a 4 x 4 matrix of plots
for i = 1:4                                       % label the plots
  xlabel(ax(4,i), labels{i});
  ylabel(ax(i,1), labels{i});
end

figure;
scatter3(data(:,1),data(:,2),data(:,3),40,data(:,4),'filled')    % draw the scatter plot
xlabel("x_1"); ylabel("x_2"); zlabel("x_3"); title('Color is x_4')

%% Possibl matrices

pts = 10000;
data = rand(pts, 4);
figure;
labels = ["x_1","x_2","x_3","x_4"];
[h,ax] = plotmatrix(data);                        % create a 4 x 4 matrix of plots
for i = 1:4                                       % label the plots
  xlabel(ax(4,i), labels{i});
  ylabel(ax(i,1), labels{i});
end

figure;
scatter3(data(:,1),data(:,2),data(:,3),40,data(:,4),'filled')    % draw the scatter plot
xlabel("x_1"); ylabel("x_2"); zlabel("x_3"); title('Color is x_4')

