
clear all; close all; clc;

%% NTKPRecursive
n = 16;
A = rand(n^2,n^2,n^2);
[B, C] = NTKP(A);
[B1, B2] = NTKP(B);
[C1, C2] = NTKP(C);

% Error
norm(A, 'fro')
norm(A - superkron(B,C), 'fro')
norm(A - superkron(superkron(B1,B2),superkron(C1,C2)), 'fro')

% Compression
prod(size(A))
prod(size(B)) * 2
prod(size(B1)) * 4


%% NTKPRank
n = 10;
A = rand(n^2,n^2,n^2);
[K, N1] = NTKPR(A, 500);

[B,sigmas] = tkpsvd(A,n * ones(1,6));
N2 = zeros(501,1); N2(1) = norm(A, 'fro');
for i=1:500
    A = A - sigmas(i) * superkron(B{1,i}, B{2,i});
    N2(i+1) = norm(A, 'fro');
end


figure; title('Multi-rank Norm fro'); hold on;
plot(N1); plot(N2); legend(['Mine', 'Yours'])

%% Identify least connected neurons to reduce the set to 16 (4x4)
D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');
rngs = [1 1508;
        1509 3017;
        3018 size(D,1)];

D1 = D(rngs(1,1):rngs(1,2),:);
D2 = D(rngs(2,1):rngs(2,2),:);
D3 = D(rngs(3,1):rngs(3,2),:);

[M1, idxs] = HAT.multicorrelations(D1, 3, 'Wang');
[M2, ~] = HAT.multicorrelations(D2, 3, 'Wang');
[M3, ~] = HAT.multicorrelations(D3, 3, 'Wang');

t = 0.9; % Disconnected near 8355
IM1 = HAT.hyperedge2IM(idxs(find(M1>t),:));
HG1 = Hypergraph('IM', IM1);
IM2 = HAT.hyperedge2IM(idxs(find(M2>t),:));
HG2 = Hypergraph('IM', IM2);
IM3 = HAT.hyperedge2IM(idxs(find(M3>t),:));
HG3 = Hypergraph('IM', IM3);

figure;
subplot(3,1,1); HG1.plot();
subplot(3,1,2); HG2.plot();
subplot(3,1,3); HG3.plot();

deg = HG1.nodeDegrees + HG2.nodeDegrees + HG2.nodeDegrees;
[v, i] = mink(deg, 5);
neurons = true(21,1); neurons(i) = false;
% Dn = D(:,neurons);

D1 = D1(:,neurons);
D2 = D2(:,neurons);
D3 = D3(:,neurons);

[M1, idxs] = HAT.multicorrelations(D1, 3, 'Wang');
[M2, ~] = HAT.multicorrelations(D2, 3, 'Wang');
[M3, ~] = HAT.multicorrelations(D3, 3, 'Wang');

t = 0.9; % Disconnected near 8355
IM1 = HAT.hyperedge2IM(idxs(find(M1>t),:));
HG1 = Hypergraph('IM', IM1); A1 = HG1.adjTensor;
IM2 = HAT.hyperedge2IM(idxs(find(M2>t),:));
HG2 = Hypergraph('IM', IM2); A2 = HG2.adjTensor;
IM3 = HAT.hyperedge2IM(idxs(find(M3>t),:));
HG3 = Hypergraph('IM', IM3); A3 = HG3.adjTensor;


[B1, C1] = NTKP(A1);
[B2, C2] = NTKP(A2);
[B3, C3] = NTKP(A3);

mean(mean(mean(B1)))
median(median(median(B1)))
figure; hold on; title('Phase 1'); histogram(B1); histogram(C1); 
figure; hold on; title('Phase 2'); histogram(B2); histogram(C2); 
figure; hold on; title('Phase 3'); histogram(B3); histogram(C3); 
