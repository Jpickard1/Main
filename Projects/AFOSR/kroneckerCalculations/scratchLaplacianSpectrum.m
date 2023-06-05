%% Laplacian Eigenspectrum
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 29, 2023


n = 7;
B = rand(n,n,n);
[b,c] = heig(B)

heig(B)

%%

n = 4;
k = 3;
HG1 = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) / 2), k);
HG2 = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) / 2), k);
HGk = Hypergraph('IM', HAT.A32IM(superkron(HG1.adjTensor, HG2.adjTensor)));

L1 = HG1.laplacianTensor;
L2 = HG2.laplacianTensor;
Lk = HGk.laplacianTensor;

[val1, vec1] = heig(L1);
[val2, vec2] = heig(L2);
[valk, veck] = heig(Lk);

%% Recreation of the first plot in Estimation of Laplacian Spectra of 
%   Direct and Strong Graph Products
clear

n = 50;
p = 0.0816;

A1 = generate_random_network('ER', n, p);
% A1 = rand(n,n); A1 = A1 + A1';
D1 = diag(sum(A1)); L1 = D1 - A1;
A2 = generate_random_network('ER', n, p);
% A2 = rand(n,n); A2 = A2 + A2';
D2 = diag(sum(A2)); L2 = D2 - A2;
Ak = kron(A1, A2); Dk = diag(sum(Ak)); Lk = Dk - Ak;
[val1, vec1] = eig(L1);
[val2, vec2] = eig(L2);
kronVecs = kron(vec1, vec2);

R = zeros(size(kronVecs,2), 1);
for i=1:size(kronVecs, 2)
    v = kronVecs(:,i);
    vout = Lk * v;
    %if sum(v) ~= 0; v = v / sum(v); end
    %if sum(vout) ~= 0; vout = vout / sum(vout); end
    r = corrcoef(v, vout);
    R(i) = r(1,2);
end

figure; plot(sort(R))
figure; histogram(R)

%% Attempt to recreate a similar experiment for hypergraphs

clear

n = 4;
k = 3;

HG1 = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) / 2), k);
HG2 = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) / 2), k);
HGk = Hypergraph('IM', HAT.A32IM(superkron(HG1.adjTensor, HG2.adjTensor)));

L1 = HG1.laplacianTensor;
L2 = HG2.laplacianTensor;
Lk = HGk.laplacianTensor;

[val1, vec1] = heig(L1);
[val2, vec2] = heig(L2);

kronVecs = kron(vec1, vec2);
kronVals = kron(val1, val2);

R = zeros(size(kronVecs,2), 1);
for i=1:size(kronVecs, 2)
    v = kronVecs(:,i);
    vout = ttvk2(tensor(Lk), kronVecs(:,i), 1);
    %if sum(v) ~= 0; v = v / sum(v); end
    %if sum(vout) ~= 0; vout = vout / sum(vout); end
    r = corrcoef(v * kronVals(i), vout);
    R(i) = r(1,2);
end

figure; plot(sort(R))
figure; histogram(R)

%% Scratch HG dynamics with regularization

clear; clc; close all;

n = 5;
k = 3;
HG = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) / 2), k);
A = HG.adjTensor;

for i=1:n
    for j=1:n
        for k=1:n
            if A(i,j,k) ~= 0 && rand(1,1) < 0.6
                A(i,j,k) = -A(i,j,k);
            end
        end
    end
end

T = 400;
X = zeros(T,n);
X(1,:) = 2 * rand(n,1) - 1;
X(1,:) = (X(1,:)); % / sum(X(1,:));
for t=2:T
    x = X(t-1,:);
    X(t,:) = (X(t-1,:)' + 0.01 * ttvk2(tensor(A), x', 1)); %reshape(A, [n, n^(k-1)]) * kron(x,x);
    % X(t,:) = X(t-1,:)' + 0.01 * reshape(A, [n, n^(k-1)]) * kron(x,x);
end

figure;
plot(X)