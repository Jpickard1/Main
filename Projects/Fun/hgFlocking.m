%% Hypergraph FLocking
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 26, 2023

% TODO: Set parameters
n = 7;
v = 4;
k = 3;
T = 10000;

% Make hypergraph
HG = HAT.uniformErdosRenyi(n,v,k);
disp(full(HG.IM));

% Get Laplacian tensor
L = HG.laplacianTensor;
L = tensor(L);

X = zeros(T,n);
x = rand(n,1);
X(1,:) = x';
for i=2:T
    v = -0.01 * ttvk(L, X(i-1,:)');
    X(i,:) = X(i-1,:) + v';
end
figure; title('State Convergence - HG'); hold on;
for i=1:n
    plot(1:T,X(:,i));
end

% Laplacian Matrix
C = HG.cliqueGraph;
L = full(diag(sum(C)) - C);

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