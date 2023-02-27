%% Hypergraph Flocking
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 26, 2023

% TODO: Set parameters
n = 10;
v = 4;
k = 3;
T = 2000;

% Make hypergraph
HG = hyperring(n,k);
% HG = HAT.uniformErdosRenyi(n,v,k);
disp(full(HG.IM));

% Get Laplacian tensor
L = HG.laplacianTensor;
% L = HG.adjTensor;
L = tensor(L);

X = zeros(T,n);
x = rand(n,1);
disp(sum(x)); disp(mean(x));
% x = ones(n,1); x(1) = 1.05; x(5) = 1.95;
X(1,:) = x';
for i=2:T
    v = -0.01 * ttvk(L, X(i-1,:)');
    X(i,:) = X(i-1,:) + v';
end
figure; title('State Convergence - HG'); hold on;
for i=1:n
    plot(1:T,X(:,i));
end

%% X is 2D cords
% TODO: Set parameters
n = 10;
v = 4;
k = 3;
T = 2000;

% Make hypergraph
HG = hyperring(n,k);
% HG = HAT.uniformErdosRenyi(n,v,k);
disp(full(HG.IM));

% Get Laplacian tensor
L = HG.laplacianTensor;
% L = HG.adjTensor;
L = tensor(L);

X = zeros(T,n,2);
x = rand(n,2);
disp(sum(x)); disp(mean(x));
% x = ones(n,1); x(1) = 1.05; x(5) = 1.95;
X(1,:,:) = x;
for i=2:T
    xx = X(i-1,:,:);
    xx = reshape(xx, 10,2 ,1);
    v = -0.01 * ttvk(L, xx);
    X(i,:) = X(i-1,:) + v';
end
figure; title('State Convergence - HG'); hold on;
for i=1:n
    plot(1:T,X(:,i));
end

%% Laplacian Matrix
C = HG.cliqueGraph;
L = full(C);
% L = full(diag(sum(C)) - C);

X = zeros(T,n);
x = rand(n,1);
X(1,:) = x';
for i=2:T
    v = -0.01 * L * X(i-1,:)';
    X(i,:) = X(i-1,:) + v';
end
figure; title('State Convergence - G'); hold on;
for i=1:n
    plot(1:T,X(:,i));
end