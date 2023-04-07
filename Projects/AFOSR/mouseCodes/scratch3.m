% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');
n = size(D, 2);

rngs = [1 1508;
        1509 3017;
        3018 size(D,1)];

% neurons = ones(15, 1); % neurons([8, 10, 13],:) = 0;
neurons = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
neurons = logical(neurons);
D1 = D(rngs(1,1):rngs(1,2), neurons);
D2 = D(rngs(2,1):rngs(2,2), neurons);
D3 = D(rngs(3,1):rngs(3,2), neurons);
[M1, idxs] = HAT.multicorrelations(D1, 3, 'Wang');
[M2, ~] = HAT.multicorrelations(D2, 3, 'Wang');
[M3, ~] = HAT.multicorrelations(D3, 3, 'Wang');
% figure; hold on; histogram(M1, 20); histogram(M2, 20); histogram(M3, 20);
t = 0.5; % Disconnected near 8355
hyperedges1 = idxs(find(M1>t),:);
IM1 = HAT.hyperedge2IM(hyperedges1);
HG1 = Hypergraph('IM', IM1);
hyperedges2 = idxs(find(M2>t),:);
IM2 = HAT.hyperedge2IM(hyperedges2);
HG2 = Hypergraph('IM', IM2);
hyperedges3 = idxs(find(M3>t),:);
IM3 = HAT.hyperedge2IM(hyperedges3);
HG3 = Hypergraph('IM', IM3);

% HOSVD decomposition
A1 = tensor(HG1.adjTensor);
S1 = hosvd(A1, 1e-12);
A2 = tensor(HG2.adjTensor);
S2 = hosvd(A2, 1e-12);
A3 = tensor(HG3.adjTensor);
S3 = hosvd(A3, 1e-12);

% Extract singular values
sv1 = [];
sv2 = [];
sv3 = [];
for i=1:size(S1.core,1)
    sv1 = [sv1 S1.core(i,i,i)];
end
for i=1:size(S2.core,1)
    sv2 = [sv2 S2.core(i,i,i)];
end
for i=1:size(S3.core,1)
    sv3 = [sv3 S3.core(i,i,i)];
end

figure; hold on;
plot(abs(sv1)); plot(abs(sv2)); plot(abs(sv3));
legend(["Feed", "Fast", "Re-Feed"]);

sv1 = [];
for i=1:size(S.core,1)
    sv1 = [sv1 S.core(i,i,i)];
end
figure; plot(sv1);
figure; plot(abs(sv1));
