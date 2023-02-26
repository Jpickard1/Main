%% Hypergraph FLocking
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 26, 2023

% TODO: Set parameters
n = 7;
v = 4;
k = 3;
T = 100;

% Make hypergraph
HG = HAT.uniformErdosRenyi(n,v,k);
disp(full(HG.IM));

% Get Laplacian tensor
L = HG.laplacianTensor;
L = tensor(L);

X = zeros(T,n);

x = rand(n,1);
X(1,:) = x;
for i=2:T
    x = ttvk(L, x);
    x = x(:);
    X(i,:) = x;
end

figure; title('State Convergence'); hold on;
for i=1:n
    plot(1:T,X(i,:));
end


