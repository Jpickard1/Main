% Auth: Joshua Pickard
%       jpic@umich.edu
% Date" April 2023

A = ones(3,3) - diag(ones(3,1));
xNames = ["x_1", "x_2", "x_3"];
yNames = ["y_1", "y_2", "y_3"];

K = kron(A,A);
kNames = kronString(xNames, yNames);

figure;
subplot(1,2,1); P = plot(graph(A)); P.NodeLabel = xNames; title('Graph');
subplot(1,2,2); P = plot(graph(K)); P.NodeLabel = kNames; title('Kronecker Graph');

x1 = sym("x_%d", [3,1]);
y1 = sym("y_%d", [3,1]);
k = kron(x1, y1);

Lx = diag(sum(A)) - A;
Lk = kron(Lx, Lx)
L = diag(sum(K)) - K







dx = Lx * x1
dy = Lx * y1
dk = Lk * k

kron(dx, dy)
