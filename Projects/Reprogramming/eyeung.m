%% This file utilizes data and methods from E. Yeung
%
% Data files were obtained from Cooper
%   - time_series_rna_fold_changes.csv
%   - time_series_rna_Z_scores.csv
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 3, 2023

%% Read in data
% Data matrices D1,D2, and D3 are replicate experiments with 15 time points
% each. The i,j-th entry the the value of RNA for gene i at time j computed
% according to the fold changes.

T = readtable("time_series_rna_fold_changes.csv");
D1 = zeros(15,size(T,2)-1);
D2 = zeros(15,size(T,2)-1);
D3 = zeros(15,size(T,2)-1);
for t=1:15
    D1(t,:) = T{3*(t-1) + 1, 2:size(T,2)};
    D2(t,:) = T{3*(t-1) + 2, 2:size(T,2)};
    D3(t,:) = T{3*(t-1) + 3, 2:size(T,2)};
end

D1 = D1'; D2 = D2'; D3 = D3';

D = [D1 D2 D3];

clearvars -except D

%% Exact DMD
A1 = exactDMD(D1);
A2 = exactDMD(D2);
A3 = exactDMD(D3);

A = A1;
X = D1(:,1);

%% Plot time scale interaction distributions
% figure; hold on; plot(sort(A1(:))); plot(sort(A2(:))); plot(sort(A3(:)));

%% Sensor Selection Problem:
%
%   dx/dt=Ax
%   y    =Cx
%
%   Maximize the output energy subject to the constraint that the output
%   matrix observes individual vertices i.e. C is spare with at most one
%   nonzero of 1 entry per row. This is expressed as:
%
%       max_C sum_{i=0}^T x0'(A^i)'C'CA^ix0 subject to CC'=I
%
%   This problem has an analytical solution from the Lagrange Dual problem.
%   In particular, the columns of C are selected to be the leading
%   eigenvectors of the observability Grammarian.
%

t = 15;
% Define the objective function
objective = @(C) -sum(arrayfun(@(i) X'*(A^i)'*C'*C*A^i*X, 1:t));

% Define the constraint function for C*C^T = I
constraint = @(C) norm(C*C' - eye(size(C,1)));

% Initial guess for C (you can provide a different initial guess)
initialC = randn(size(X,1), size(X,1));

% Set up optimization options
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'MaxIterations', 1000);

% Solve the optimization problem
[C_opt, fval] = fmincon(objective, initialC, [], [], [], [], [], [], constraint, options);

% Display the optimized C
disp('Optimized C:');
disp(C_opt);

% Display the maximum value of the objective function
disp(['Maximum value of the objective function: ', num2str(-fval)]);


%% Observability Grammarian



%% Kalman Filter

C1 = zeros(1,size(A1,1)); C1(1) = 1;

Ts = -1;
sys = ss(A1,[],C1,[],Ts,'OutputName','y');  % Plant dynamics and additive input noise w
[kalmf,L,~,Mx,Z] = kalman(sys,[],[]);


%% Sensor Selection


%% Apply DMD

O1 = DMD(D1, [], 0.9)
O1.DMD
