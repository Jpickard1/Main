%% Hypergraph Agreement Protocol
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 31, 2023

clear

n = 5;
k = 3;
HG = HAT.uniformErdosRenyi(n,round(0.5 * nchoosek(n,k)),k);

L = HG.laplacianTensor();

T = 100; n1=n;
s = 2;
s1 = 0.01;
s2 = 0.01;
t1 = [0:(T-1)] * s1;
t2 = [0:(T-1)] * s2;

X1 = zeros(T, n1);  X1(1,:) = s * rand(n1,1) - (s/2);
for t=2:T
    X1(t,:) = X1(t-1,:)' + s1 * ttvk(tensor(L), X1(t-1,:)');
end
figure; plot(t1, X1);
figure; HG.plot()
zeig(L)

%%
X2 = zeros(T, n1^2);  X2(1,:) = kron(X1(1,:),X1(1,:)); %s * rand(n1^2,1) - (s/2);
for t=2:T
    X2(t,:) = X2(t-1,:)' + s2 * ttvk(tensor(B), X2(t-1,:)');
end

figure;
subplot(1,2,1); plot(t1, X1); ylabel('State'); xlabel('Time'); title('')
subplot(1,2,2); plot(t2, X2); ylabel('State'); xlabel('Time'); title('Kronecker System')