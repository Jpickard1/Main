%% 03/06/2023 - Are Linear Time Invariant Systems observable up to a scaler output?
%
%   It would seem so.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 6, 2023

%% Scaled observability test
%   In this test, we construct a the system (1) and the corresponding
%   Luendberg full state estimator.
%
%       (System 1)  dx/dt = Ax, y=Cx
%
%   When y=Cx, the error of the estimator decays to 0 over time. However,
%   we define lambda as a random constant, and take output y=lambda*C*x.
%   Then, we define the error of the estimator with equation 2
%
%       (Equation 2)    epsion = lambda * X - Xhat
%
%   Numerically, we see that the error function in Equation 2 goes to 0
%   over time.

clear; clc; close all;

%TODO: SET PARAMETERS
RNG = 100; % Initial random conditions
lambda = rand(); % Set scaling
T = 10000;
plt = true;

% SYSTEM PARAMETER
A = zeros(2,2);
A(1,1) = 1.03;      % Set eigenvalues of the system
A(2,2) = 0.97;
P = rand(2,2);
A = 0.001 * P * A * inv(P); % State transition matrix
C = [1 0];          % Observability Matrix

disp(eig(A))
disp(A)
if rank(obsv(A,C)) < 2
    disp('NOT OBSERVABLE');
end

% SYSTEM MODEL
X = zeros(T,2);
X(1,:) = RNG * rand(2,1);
for t=2:T
    X(t,:) = X(t-1,:)' + A * X(t-1,:)';
end
systemOutput = lambda * C * X';

% FULL STATE OBSERVER
% Kf = lqe(A,zeros(size(A)),C,zeros(size(A)),ones([1 1]),zeros([2 1]));
Kf = (lqr(A',C',eye(size(A)),eye([1 1])))';
Xhat = zeros(T,2);
Xhat(1,:) = RNG * rand(2,1);
systemOutputHat = zeros(size(systemOutput));
for t=2:T
    Xhat(t,:) = Xhat(t-1,:)' + A * Xhat(t-1,:)' + Kf * (systemOutput(t-1) - systemOutputHat(t-1));
    systemOutputHat(t) = C * Xhat(t,:)';
end

epsilon = lambda * X - Xhat;

if plt
    figure; hold on; title('System and Estimator Variables');
    plot(X(:,1));
    plot(X(:,2));
    plot(Xhat(:,1));
    plot(Xhat(:,2));
    legend(["x1", "x2", "x1hat", "x2hat"]);
    xlabel('Time'); ylabel('State');

    figure; hold on; title("Norm of Epsion (lambda=" + string(lambda) + ")");
    plot(sqrt(sum(epsilon.^2,2)));
    xlabel('Time'); ylabel('Error');
end

%% Kronecker Controllability Test
%   In this test we consturct a network A = kron(A1,A2). We construct
%   matrices B1 and B2 so (Ai, Bi) are controllable. Then, we test if 
%   (A, kron(B1, B2)) is controllable. This is to numerically validate
%   Theorem 6 in Kronecker Product of Networked Systems and their
%   Approximates.
clear; clc; close all;

% TODO: Set params
x = 1;
while x<100
x = x + 1;
n1 = 4; % Size of A1
n2 = 4; % Size of A2
p1 = 1; % B1 is n1 x p1
p2 = 1; % B2 is n2 x p2

A1 = erdos_renyi_network(n1, round(0.5*nchoosek(n1,2)));
A2 = erdos_renyi_network(n2, round(0.5*nchoosek(n1,2)));
A1 = double(A1); A2 = double(A2);
B1 = zeros(n1, p1); vxc1 = randi([1, n1], p1, 1);
for i=1:length(vxc1); B1(vxc1(i),i) = 1; end
B2 = zeros(n2, p2); vxc2 = randi([1, n2], p2, 1);
for i=1:length(vxc2); B2(vxc2(i),i) = 1; end
r1 = rank(ctrb(A1,B1));
r2 = rank(ctrb(A2,B2));
if r1 ~= n1 || r2 ~= n2
    disp('The small systems are not controllable');
    continue;
end

A = kron(A1, A2); B = kron(B1, B2);
r =  rank(ctrb(A,B));

% If (A, C) is not controllable/observable
if r ~= (n1 * n2)
    disp("Not Observable")
    [v, d] = eig(A);
    isD = isdiag(inv(v) * A * v);
    
    [~,~,W1] = eig(A1);
    [~,~,W2] = eig(A2);
    s1 = W1' * B1;
    s2 = W2' * B2;

    if ~isD
        disp("A is not diagnolizable");
    elseif (sum(find(s1 == 0)) > 0) || (sum(find(s2 == 0)) > 0)
        disp("Condition 2");
    elseif isD 
        disp("THM is Wrong!");
    end
end
end

%% Kronecker Product of Networked Systems and their Approximates Example 1
% 
%   Here I verify Example 7's comments on the Kronecker controllability of
%   the graphs in Example 1.

clear; close all; clc

A1 = [1 1 0;
      0 0 1;
      1 0 0];
A2 = [1 1 0;
      1 1 1;
      0 1 1];
B1 = [1; 0; 0]; C1 = B1';
B2 = [1; 1; 0]; C2 = B2';

c1 = rank(ctrb(A1,B1));
o1 = rank(obsv(A1,C1));
c2 = rank(ctrb(A2,B2));
o2 = rank(obsv(A2,C2));

A = kron(A1, A2);
B = kron(B1, B2);
C = kron(C1, C2);

c = rank(ctrb(A,B))
o = rank(obsv(A,C))

%% Full State Estimator for Kronecker Observed Hypergraphs Try - 1
%
%   Here I try and design a Leunberg observer for observable Kronecker
%   graphs

clear; close all; clc

% TODO: Set parameters
T = 1000; % Time
RNG = 10;
n1 = 4; % Size of A1
n2 = 4; % Size of A2
p1 = 2; % B1 is n1 x p1
p2 = 2; % B2 is n2 x p2

% 1. CONSTRUCT SYSTEM
% Select kronecker observable hypergraphs
O = 1; 
while O > 0
    O = O + 1;
    % Construct networks
    A1 = erdos_renyi_network(n1, round(0.5*nchoosek(n1,2))); A1 = real(A1);
    A2 = erdos_renyi_network(n2, round(0.5*nchoosek(n1,2))); A2 = real(A2);
    
    % Construct observation matrices
    C1 = zeros(p1, n1); vxc1 = randi([1, n1], p1, 1);
    for i=1:length(vxc1); C1(i,vxc1(i)) = 1; end
    C2 = zeros(p2, n2); vxc2 = randi([1, n2], p2, 1);
    for i=1:length(vxc2); C2(i,vxc2(i)) = 1; end

    A = kron(A1, A2); C = kron(C1, C2);
    if rank(obsv(A,C)) == n1 * n2
        disp(O); O = 0;
    end
end

%% 2. SYSTEM MODEL
X = zeros(T,n1*n2);
X0a1 = RNG * rand(n1,1);
X0a2 = RNG * rand(n2,1);
X(1,:) = kron(X0a1, X0a2);
% X(1,:) = RNG * rand(n1*n2,1);
for t=2:T
    X(t,:) = X(t-1,:)' + A * X(t-1,:)';
end
systemOutput = C * X';

%% 3. FACTOR SYSTEM OUTPUT TO CONSTRUCT OUTPUTS FOR A1 AND A2
systemOutputReshaped = reshape(systemOutput, 2,2, T);
pseudoA1output1 = zeros(p1, T); % Each subsystem gets multiple outputs because SVD has multiple singular vectors
pseudoA2output1 = zeros(p2, T);
pseudoA1output2 = zeros(p1, T);
pseudoA2output2 = zeros(p2, T);
for t=1:T
    sot = systemOutputReshaped(:,:,t);
    [U,S,V] = svd(sot);
    pseudoA1output1(:,t) = U(:,1);
    pseudoA2output1(:,t) = V(:,1);
    pseudoA1output2(:,t) = U(:,2);
    pseudoA2output2(:,t) = V(:,2);
    % U(:,1) * (S(1,1)) * V(:,1)'; % + U(:,2) * (S(2,2)) * V(:,2)'
end

% Full state estimators for A1 and A2
X1hat1 = fse(A1, C1, pseudoA1output1, RNG);
X2hat1 = fse(A2, C2, pseudoA2output1, RNG);
X1hat2 = fse(A1, C1, pseudoA1output2, RNG);
X2hat2 = fse(A2, C2, pseudoA2output2, RNG);

Xhat = 

%% 4. Plot all 16 trajectories







%%
c = kroneckerObservable(A1,A2,C1,C2)

A = kron(A1, A2); C = kron(C1, C2); O = obsv(A, C);
rank(O)

%% Kronecker Observability & Controllability Test
%   In this test we consturct a network A = kron(A1,A2). We construct
%   matrices C1 and C2 so (Ai, Ci) are observable. Then, we test if 
%   (A, kron(C1, C2)) is observable.
clear; clc; close all;

% TODO: Set params
n1 = 4; % Size of A1
n2 = 4; % Size of A2
p1 = 1; % C1 is n1 x p1
p2 = 1; % C2 is n2 x p2

A1 = rand(n1,n1); A1 = A1 + A1';
A2 = rand(n2,n2); A2 = A2 + A2';
C1 = zeros(p1, n1); vxc1 = randi([1, n1], p1, 1);
for i=1:length(vxc1); C1(i,vxc1(i)) = 1; end
C2 = zeros(p2, n2); vxc2 = randi([1, n2], p2, 1);
for i=1:length(vxc2); C2(i,vxc2(i)) = 1; end

A = kron(A1, A2); C = kron(C1, C2);
%{
r1 = rank(obsv(A1,C1));
r2 = rank(obsv(A2,C2));
r = rank(obsv(A,C));
%}
r1 = rank(ctrb(A1,C1'));
r2 = rank(ctrb(A2,C2'));
r =  rank(ctrb(A,C'));

if rank(C1) < p1
    disp("===")
elseif rank(C2) < p2
    disp("===")
elseif r1 < n1
    disp("===")
elseif r2 < n2
    disp("===")
end

disp(r1);
disp(r2);
disp(r);

[W, E] = eig(A1);
disp(W * C1')

[W, E] = eig(A2);
disp(W * C2')

%%
% SYSTEM PARAMETER
A1 = rand(n1, n1);
A1 = zeros(2,2);
A1(1,1) = 1.03;      % Set eigenvalues of the system
A1(2,2) = 0.97;
P = rand(2,2);
A1 = 0.001 * P * A1 * inv(P); % State transition matrix
C1 = [1 0];          % Observability Matrix
if rank(obsv(A1,C1)) < 2
    disp('NOT OBSERVABLE');
end
A2 = zeros(2,2);
A2(1,1) = 1.03;      % Set eigenvalues of the system
A2(2,2) = 0.97;
P = rand(2,2);
A2 = 0.001 * P * A2 * inv(P); % State transition matrix
C2 = [1 0];          % Observability Matrix
if rank(obsv(A2,C2)) < 2
    disp('NOT OBSERVABLE');
end
A = kron(A1, A2);
C = kron(C1, C2);
if rank(obsv(A,C)) < 4
    disp('NOT OBSERVABLE');
end

%%TODO: SET PARAMETERS
RNG = 100; % Initial random conditions
lambda = 1;%rand(); % Set scaling
T = 10000;
plt = true;

% SYSTEM PARAMETER
A1 = zeros(2,2);
A1(1,1) = 1.03;      % Set eigenvalues of the system
A1(2,2) = 0.97;
P = rand(2,2);
A1 = 0.001 * P * A1 * inv(P); % State transition matrix
C1 = [1 0];          % Observability Matrix
if rank(obsv(A1,C1)) < 2
    disp('NOT OBSERVABLE');
end
A2 = zeros(2,2);
A2(1,1) = 1.03;      % Set eigenvalues of the system
A2(2,2) = 0.97;
P = rand(2,2);
A2 = 0.001 * P * A2 * inv(P); % State transition matrix
C2 = [1 0];          % Observability Matrix
if rank(obsv(A2,C2)) < 2
    disp('NOT OBSERVABLE');
end
A = kron(A1, A2);
C = kron(C1, C2);
if rank(obsv(A,C)) < 4
    disp('NOT OBSERVABLE');
end

%% SYSTEM MODEL
X = zeros(T,2);
X(1,:) = RNG * rand(2,1);
for t=2:T
    X(t,:) = X(t-1,:)' + A * X(t-1,:)';
end
systemOutput = lambda * C * X';

% FULL STATE OBSERVER
% Kf = lqe(A,zeros(size(A)),C,zeros(size(A)),ones([1 1]),zeros([2 1]));
Kf = (lqr(A',C',eye(size(A)),eye([1 1])))'
Xhat = zeros(T,2);
Xhat(1,:) = RNG * rand(2,1);
systemOutputHat = zeros(size(systemOutput));
for t=2:T
    Xhat(t,:) = Xhat(t-1,:)' + A * Xhat(t-1,:)' + Kf * (systemOutput(t-1) - systemOutputHat(t-1));
    systemOutputHat(t) = C * Xhat(t,:)';
end

epsilon = lambda * X - Xhat;

if plt
    figure; hold on; title('System and Estimator Variables');
    plot(X(:,1));
    plot(X(:,2));
    plot(Xhat(:,1));
    plot(Xhat(:,2));
    legend(["x1", "x2", "x1hat", "x2hat"]);
    xlabel('Time'); ylabel('State');

    figure; hold on; title('Norm of Epsion');
    plot(sqrt(sum(epsilon.^2,2)));
    xlabel('Time'); ylabel('Error');
end