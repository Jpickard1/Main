
clear; close all; clc


%% Movement model
a = -0.5
b = 0.1
c = 0.2

M = zeros(2,2,2);

M(1,1,1) = a;
M(2,1,1) = -a;
M(1,2,1) = b;
M(2,1,2) = -b;
M(1,2,2) = c;
M(2,2,2) = -c;

T = 400;
X = zeros(T,2);
X(1,:) = rand(2,1); X(1,:) = X(1,:) / sum(X(1,:));
% X(1,:) = [0.8410    0.1590];
for t=2:T
    x = X(t-1,:);
    X(t,:) = X(t-1,:)' + 0.01 * reshape(M, [2 4]) * kron(x,x)';
end

figure; hold on;
plot(X(:,1));
plot(X(:,2));
xlabel('Time');
ylabel('Population Size');
title('Inter Population Dynamics');
legend(["Population 1", "Population 2"]);

% Good initial conditions
%     0.5999    0.4001
%     0.5783    0.4217
%     0.8410    0.1590

%% Lotka-Voltera Model

alpha = 0.4;
beta = 0.01;
gamma = 0.02;
delta = 0.01;

T = 40;
X = zeros(T,2);
X(1,:) = [20 20];
for t=2:T
    X(t,1) = X(t-1,1) + 0.01 * (alpha * X(t-1,1) - beta * X(t-1,1) * X(t-1,2));
    X(t,2) = X(t-1,2) + 0.01 * (delta * X(t-1,1) * X(t-1,2) - gamma * X(t-1,2));
end

figure;
plot(X)
hold on;
legend(["Prey 1", "Predator 2"]);

%% Quadratic Lotka-Voltera Model

alpha = 0.5;
beta  = 0.02;
gamma = 0.02;
delta = 0.03;

T = 6000;
X = zeros(T,2);
X(1,:) = [0.5 0.5];
for t=2:T
    X(t,1) = X(t-1,1) + 0.01 * (alpha * X(t-1,1) - beta * X(t-1,1) * X(t-1,2));
    X(t,2) = X(t-1,2) + 0.01 * (beta * X(t-1,1) * X(t-1,2) - gamma * X(t-1,2) * X(t-1,2));
    % X(t,2) = X(t-1,2) + 0.01 * (delta * X(t-1,1) * X(t-1,2) - gamma * X(t-1,2) * X(t-1,2));
    if sum(X(t,:)) > 200; break; end
end

figure;
plot(X)
hold on;
legend(["Prey", "Predator"]);


%%
A(1,1,1) = alpha; % a;
A(2,1,1) = alpha; % -a;
A(1,2,1) = beta;  %b;
A(2,1,2) = beta;  %-b;
A(1,2,2) = c;
A(2,2,2) = -c;

%% MATLAB Lotka-Voltera Model

t0 = 0;
tfinal = 15;
y0 = [20; 20];   
[t,y] = ode23(@lotka,[t0 tfinal],y0);


plot(t,y)
title('Predator/Prey Populations Over Time')
xlabel('t')
ylabel('Population')
legend('Prey','Predators','Location','North')

plot(y(:,1),y(:,2))
title('Phase Plane Plot')
xlabel('Prey Population')
ylabel('Predator Population')


%% Movement model with ode45

tspan = [0 20];
y0 = [0.8410    0.1590];
sol = ode45(@mmodel ,tspan,y0)

figure; hold on;
plot(sol.x, sol.y(1,:))
plot(sol.x, sol.y(2,:))

%{
function dxdt = mmodel(t,x)
    a = 0.5;
    b = 0.3;
    c = 0.2;
    dxdt = zeros(2,1);
    dxdt(1) =   a * x(1)^2 + b * x(1) * x(2) + c * x(2)^2;
    dxdt(2) = - a * x(1)^2 - b * x(1) * x(2) - c * x(2)^2;
end
%}

%% oscilating hypergraph dynamics

clear; clc; close all;

n = 5;
k = 4;
HG = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) / 2), k);
A = HG.adjTensor;

T = 400;
X = zeros(T,n);
X(1,:) = rand(n,1);
X(1,:) = X(1,:) / sum(X(1,:));
for t=2:T
    x = X(t-1,:);
    X(t,:) = X(t-1,:)' + 0.01 * reshape(A, [n, n^(k-1)]) * kron(x,x);
end

figure;
plot(X)

%% Laplacian heig vectors

A = rand(3,3,3);
[b,c] = heig(A)

%%

sum(sum(sum(sum(A))))

for i=1:size(A,1)
    for j=1:size(A,1)
        for k=1:size(A,1)
            for l=1:size(A,1)
                if k < i || l < j
                    A(i,j,k,l) = -A(i,j,k,l);
                end
            end
        end
    end
end
sum(sum(sum(sum(A))))




