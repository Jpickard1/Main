%% Scratch

%% September 10, 2023
%
%   Here we perform the following:
%     SYNTHETIC SYSTEM
%       1. construct a directed, weighted adjacency matrix A
%       2. simulate time series data
%     DMD Identification
%       1. Perform exact DMD based on the time series data to identify A
%     STRUCTURE-ID
%       1. Ap = unweighted, undirected adjacency matrix of A
%       2. Approximate the parameters of A using the sparse structure and
%       the time series data
%     COMPARISON
%       1. On a new set of time series data produced from A, I will compare
%       the predictive ability of Ad (DMD) and As (STRUCTURE ID)

%%     SYNTHETIC SYSTEM
%       1. construct a directed, weighted adjacency matrix A
n = 10;     % number of states
d = 0.25;   % density of system
A = sprand(n,n,d);
A = full(A);
%       2. simulate time series data
f = @(t,x) A * x;
[t, X] = ode45(f, [0, 100], rand(n,1));

figure;
plot(X)

%%     DMD Identification
%       1. Perform exact DMD based on the time series data to identify A
Ad = exactDMD(X');

%%     STRUCTURE-ID
%       1. Ap = unweighted, undirected adjacency matrix of A
%       2. Approximate the parameters of A using the sparse structure and
%       the time series data
Ap = double(A ~= 0);

%% GTP-4.0 Example

% A0 = A;
dXdt = X(1:size(X,1)-1,:) - X(2:size(X,1),:);
A = dXdt;
B = X(1:size(X,1)-1,:);
% A = [ ... ];  % Known matrix
% B = [ ... ];  % Known matrix
[r,c] = find(Ap == 0);
I = [c r];
n = size(B,2);  % The size of C should be compatible with B

% Formulate the Quadratic Programming problem
H = 2 * (B' * B);
f = -B' * A; 
% We flatten matrix C to a vector. Thus, we need to flatten H and f to match. 

H = kron(H,eye(n));
f = f(:);

% Setup the equality constraints 'C_ij = 0'
Aeq = zeros(size(I,1), n*n);
beq = zeros(size(I,1), 1);
for k = 1:size(I,1)
    Aeq(k, sub2ind([n, n], I(k,2), I(k,1))) = 1;
end

% Solve the QP problem
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
C = quadprog(H, f, [], [], Aeq, beq, [], [], [], options);  

% Reshape the solution vector back to a matrix form
C = reshape(C,[n,n]);
