% Trajectory approxmation of Kronekcer hypergraphs
%
% X1(1,:) = [   -0.4846   -0.4570   -0.3310    0.1491];
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 7, 2023


% Plotting commands
rPlts = 1;
cPlts = 3;

% Make hypergraphs
k = 3;
n1 = 4;
n2 = n1;
T = 100;

HG1 = getToyHG(n1,k,'hyperchain');
A1 = HG1.adjTensor; A = superkron(A1, A1);
A1 = tensor(A1); A = tensor(A);

X1 = zeros(T, n1);  X1(1,:) = rand(n1,1) - 0.5;
for t=2:T
    X1(t,:) = 0.04 * ttvk(A1, X1(t-1,:)') + X1(t-1,:)';
end
figure; subplot(rPlts,cPlts,1); hold on; plot(X1); ylabel("Hyperchain"); title("Hypergraph")
X = zeros(T, n1*n1); X(1,:) = kron(X1(1,:), X1(1,:));
for t=2:T
    X(t,:) = 0.04 * ttvk(A, X(t-1,:)') + X(t-1,:)';
end
subplot(rPlts,cPlts,2); hold on; plot(X); title("Kronecker Hypergraph")
D = 0.01 * tensor(rand(size(A)));
Xn = zeros(T, n1*n1); Xn(1,:) = kron(X1(1,:), X1(1,:));
E = zeros(T, n1*n1);
for t=2:T
    Xn(t,:) = 0.04 * ttvk(A+D, Xn(t-1,:)') + Xn(t-1,:)';
    E(t,:) = 0.04 * ttvk(D, Xn(t-1,:)');
end
N = X - Xn; N = sqrt(sum(N.^2,2)); E = sqrt(sum(E.^2,2));
subplot(rPlts,cPlts,3); hold on; plot(N); plot(E); legend(["True", "Predicted"]); title("Error")

%% 
HG1 = getToyHG(n1,k,'hyperstar');
A1 = HG1.adjTensor; A = superkron(A1, A1);
A1 = tensor(A1); A = tensor(A);
X1 = zeros(T, n1);  X1(1,:) = rand(n1,1) - 0.5;
for t=2:T
    X1(t,:) = 0.04 * ttvk(A1, X1(t-1,:)') + X1(t-1,:)';
end
subplot(rPlts,cPlts,3); hold on; plot(X1); ylabel("Hyperstar");
X = zeros(T, n1*n1); X(1,:) = kron(X1(1,:), X1(1,:));
for t=2:T
    X(t,:) = 0.04 * ttvk(A, X(t-1,:)') + X(t-1,:)';
end
subplot(rPlts,cPlts,4); hold on; plot(X);

HG1 = getToyHG(n1,k,'hyperstar');
A1 = HG1.adjTensor; A = superkron(A1, A1);
A1 = tensor(A1); A = tensor(A);
X1 = zeros(T, n1);  X1(1,:) = rand(n1,1) - 0.5;
for t=2:T
    X1(t,:) = 0.04 * ttvk(A1, X1(t-1,:)') + X1(t-1,:)';
end
subplot(rPlts,cPlts,5); hold on; plot(X1); ylabel("Hyperring");
X = zeros(T, n1*n1); X(1,:) = kron(X1(1,:), X1(1,:));
for t=2:T
    X(t,:) = 0.04 * ttvk(A, X(t-1,:)') + X(t-1,:)';
end
subplot(rPlts,cPlts,6); hold on; plot(X);

