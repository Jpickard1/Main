%% Spatial Graphs

A = uniform_spatial_graph(5000, 10, 0.2, 2, 'euclidean', true)

%% Linear Algebra check

th = 0:pi/50:2*pi;
x = 1 * cos(th);
y = 1 * sin(th);
unit_circle = [x' y'];

A = [1 2; 0 3];

P = A * unit_circle';

figure;
hold on;
plot(unit_circle(:,1), unit_circle(:,2))
plot(P(1,:), P(2,:));

[vecs,vals] = eig(A);

plot([0 1], [0 0])
plot([0 1], [0 1])


