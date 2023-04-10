% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 6, 2023

%% Testing NTKP (2)
%   X axis = dimensions of B and C
%   Y axis = error

itrs = 5;
D = zeros(20, itrs);
for n1=2:20
    for itr=1:itrs
        A1 = rand(n1, n1, n1); A2 = rand(n1, n1, n1);
        A = superkron(A1, A2);
        [B, C] = NTKP(A);
        eps = A - superkron(B, C);
        D(n1, itr) = sum(sum(sum(eps))) / sum(sum(sum(A)));
        disp(itr)
    end
    disp(n1);
end

figure; plot(mean(D,2))

%% Testing NTKP (1)

n1 = 8; k = 3;
A1 = rand(n1, n1, n1); A2 = rand(n1, n1, n1);
A = superkron(A1, A2);

[B, C] = NTKP(A);

eps = A - superkron(B, C);
sum(sum(sum(eps)))
sum(sum(sum(A)))

%% Algorithm for Kronecker Tensor Factorization
%   NOTES:
%       - currently fixed for k=3

clear; clc;
% Input: a tensor A
%{
    n1 = 8; k = 3;
    A1 = rand(n1, n1, n1); A2 = rand(n1, n1, n1);
    A = superkron(A1, A2);
%}

% Flatten tensor
Amat = reshape(A, [size(A,1), size(A,2) * size(A,3)]);

% Generate column permutations (this is the part that needs more intuition)
n = 1:n1;
p = [];
for i=1:n1
    p = [p (k*n1*(i-1) + n)];
end
rp = [];
for i=1:n1
    rp = [rp (n1*(i-1) + p)];
end
Permutations = [];
for i=1:n1
    Permutations = [Permutations ((n1^k)*(i-1) + rp)];
end
AmatP = Amat(:, Permutations);
% Nearest Kronecker Product
[Bv, Cv] = nearestKroneckerProduct(AmatP, [n1 n1^2], [n1 n1^2]);
% Output tensors
Bt = reshape(Bv, [n1, n1, n1]);
Ct = reshape(Cv, [n1, n1, n1]);

sum(sum(AmatP - kron(Bv, Cv)))

BtCt = superkron(Bt, Ct);
sum(sum(sum(A - BtCt)))
sum(sum(sum(A)))
% A = A- BtCt;
% BC = kron(B, C)


%% Tensor kronecker factors
%   It works in this section
clear; clc;

% Construct adjacency matrices
A1 = rand(3,3,3); A2 = rand(3,3,3);
A = superkron(A1, A2);

% Unfold adjacency tensors
A1mat = reshape(A1, 3, 9);
A2mat = reshape(A2, 3, 9);
Amat = reshape(A, 9, 81);

AmatK = kron(A1mat, A2mat);

rows = [1:3 10:12 19:21 4:6 13:15 22:24 7:9 16:18 25:27];
rows = [rows (27 + rows) (54 + rows)];
AM = Amat(:,rows);

% Verify if permutations are correct
E = zeros(9,81);
for i=1:9
    for j=1:81
        if ~isequal(AM(i,j), AmatK(i,j))
            E(i,j) = 1;
        end
    end
end
disp(sum(sum(E)))

[B, C] = nearestKroneckerProduct(AM, [3 9], [3 9])
BC = kron(B, C)

BC - AmatK

T = reshape(BC, [9 9 9]);

A - T

AM - BC

%% Tensor kronecker factors
clear; clc;

% Construct adjacency matrices
A1 = sym("x", [3, 3, 3]);
A2 = sym("y", [3, 3, 3]);
A = superkron(A1, A2);

% Unfold adjacency tensors
A1mat = reshape(A1, 3, 9);
A2mat = reshape(A2, 3, 9);
Amat = reshape(A, 9, 81);

AmatK = kron(A1mat, A2mat);

rows = [1:3 10:12 19:21 4:6 13:15 22:24 7:9 16:18 25:27];
rows = [rows (27 + rows) (54 + rows)];
AM = Amat(:,rows);

%S1 = shuffleMatrix(9, 9);
%SX = (Amat * S1 * S1 * S1);
%disp(SX(1,:))

%S1 = shuffleMatrix(9, 9);
%SX = (S1 * Amat')';
%disp(SX(1,:))
%disp(Amat(1,:));

% Verify if permutations are correct
E = zeros(9,81);
for i=1:9
    for j=1:81
        if ~isequal(AM(i,j), AmatK(i,j))
            E(i,j) = 1;
            % disp(string(i) + ":" + string(j) + " - " + string(Amat(i,j)) + " -- " + string(AmatK(i,j)));
        end
    end
end
disp(E)

%% Shuffle matrices

clear
X = sym('x', [9, 9]);
S = shuffleMatrix(3,3);
SX = S' * X;


XS

%% Kronecker factor 2
%   Problem: Given adjacency matrices with A1 and A2 and H with A =
%   kron(A1, A2), what is the relationship between the matriceized
%   adjacency tensors?

%{
Amat(1,:)
[x1_1_1   x1_1_1   x1_1_1
 x1_2_1   x1_2_1   x1_2_1
 x1_3_1   x1_3_1   x1_3_1
 x1_1_1   x1_1_1   x1_1_1
 x1_2_1   x1_2_1   x1_2_1
 x1_3_1   x1_3_1   x1_3_1
 x1_1_1   x1_1_1   x1_1_1
 x1_2_1   x1_2_1   x1_2_1
 x1_3_1   x1_3_1   x1_3_1

 x1_1_2   x1_1_2   x1_1_2
 x1_2_2   x1_2_2   x1_2_2
 x1_3_2   x1_3_2   x1_3_2
 x1_1_2   x1_1_2   x1_1_2
 x1_2_2   x1_2_2   x1_2_2
 x1_3_2   x1_3_2   x1_3_2
 x1_1_2   x1_1_2   x1_1_2
 x1_2_2   x1_2_2   x1_2_2
 x1_3_2   x1_3_2   x1_3_2

 x1_1_3   x1_1_3   x1_1_3
 x1_2_3   x1_2_3   x1_2_3
 x1_3_3   x1_3_3   x1_3_3
 x1_1_3   x1_1_3   x1_1_3
 x1_2_3   x1_2_3   x1_2_3
 x1_3_3   x1_3_3   x1_3_3
 x1_1_3   x1_1_3   x1_1_3
 x1_2_3   x1_2_3   x1_2_3
 x1_3_3   x1_3_3   x1_3_3

AmatK(1,:)
 x1_1_1   x1_1_1   x1_1_1
 x1_1_1   x1_1_1   x1_1_1
 x1_1_1   x1_1_1   x1_1_1
 x1_2_1   x1_2_1   x1_2_1
 x1_2_1   x1_2_1   x1_2_1
 x1_2_1   x1_2_1   x1_2_1
 x1_3_1   x1_3_1   x1_3_1
 x1_3_1   x1_3_1   x1_3_1
 x1_3_1   x1_3_1   x1_3_1

 x1_1_2   x1_1_2   x1_1_2
 x1_1_2   x1_1_2   x1_1_2
 x1_1_2   x1_1_2   x1_1_2
 x1_2_2   x1_2_2   x1_2_2
 x1_2_2   x1_2_2   x1_2_2
 x1_2_2   x1_2_2   x1_2_2
 x1_3_2   x1_3_2   x1_3_2
 x1_3_2   x1_3_2   x1_3_2
 x1_3_2   x1_3_2   x1_3_2

 x1_1_3   x1_1_3   x1_1_3
 x1_1_3   x1_1_3   x1_1_3
 x1_1_3   x1_1_3   x1_1_3
 x1_2_3   x1_2_3   x1_2_3
 x1_2_3   x1_2_3   x1_2_3
 x1_2_3   x1_2_3   x1_2_3
 x1_3_3   x1_3_3   x1_3_3
 x1_3_3   x1_3_3   x1_3_3
 x1_3_3   x1_3_3   x1_3_3
%}

clear; clc;

% Construct adjacency matrices
A1 = sym("x", [3, 3, 3]);
A2 = sym("y", [3, 3, 3]);
% A1 = rand(3,3,3);
% A2 = rand(3,3,3);
A = superkron(A1, A2);

% Unfold adjacency tensors
A1mat = reshape(A1, 3, 9);
A2mat = reshape(A2, 3, 9);
Amat = reshape(A, 9, 81);

AmatK = kron(A1mat, A2mat);

S = shuffleMatrix(3, 3);
AmatS = (Amat' * S)';
disp(AmatS(1,1:9))
disp(AmatK(1,1:9))


E = zeros(9,81);
for i=1:9
    for j=1:81
        if ~isequal(AmatS(i,j), AmatK(i,j))
            E(i,j) = 1;
            % disp(string(i) + ":" + string(j) + " - " + string(Amat(i,j)) + " -- " + string(AmatK(i,j)));
        end
    end
end
disp(E)


[B, C] = nearestKroneckerProduct(Amat, [3 9], [3 9])

% Create shuffle matrix
S = shuffleMatrix(3, 3);
R = kron(A1mat * S, A2mat * S);

sum(sum((R - Amat) > 1e-4))

for i=1:27
    disp(AmatK(1,3*(i-1)+1:3*(i-1)+3));
end

%% Kron factor
clear
HG = getToyHG(9,3,'hyperstar');
A = HG.adjTensor;
Amat = reshape(A, 9, 81);

[Bmat, Cmat] = nearestKroneckerProduct(Amat, [3 9], [3 9]);

B = reshape(Bmat, [3 3 3]);
C = reshape(Cmat, [3 3 3]);

Ap = superkron(B, C);

%% 
clear; clc;
HG1 = getToyHG(3,3,'hyperstar'); A1 = HG1.adjTensor;
HG2 = getToyHG(3,3,'hyperstar'); A2 = HG2.adjTensor;

A1mat = reshape(A1, 3, 9);
A2mat = reshape(A2, 3, 9);

S = shuffleMatrix(3, 3);

A1s = S * A1mat';
A2s = S * A2mat';

Amat = kron(A1s, A2s);
A = reshape(Amat, 9, 9, 9);


Amat = kron(A1mat, A2mat);
A = reshape(Amat, 9, 9, 9);


% [U, S, Vt] = svd(Amat);

% [O, J] = HGObsvSym(HG)

%% Kron factors
D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');
n = size(D, 2);

rngs = [1 1508;
        1509 3017;
        3018 size(D,1)];

% neurons = ones(15, 1); % neurons([8, 10, 13],:) = 0;
neurons = ones(16,1);
% neurons = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
neurons = logical(neurons);
D1 = D(rngs(1,1):rngs(1,2), neurons);
[M1, idxs] = HAT.multicorrelations(D1, 3, 'Wang');
t = 0.5; % Disconnected near 8355
hyperedges1 = idxs(find(M1>t),:);
IM1 = HAT.hyperedge2IM(hyperedges1);
HG1 = Hypergraph('IM', IM1);

% HOSVD decomposition
A1 = tensor(HG1.adjTensor);
S1 = hosvd(A1, 1e-12);

(S1.u{1} - S1.u{2}) < 1e-8

u1 = S1.u{1}(:,1)
U1 = reshape(u1,[4,4])

%% HOSVD Scree plot
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


%%
sv1 = [];
for i=1:size(S.core,1)
    sv1 = [sv1 S.core(i,i,i)];
end
figure; plot(sv1);
figure; plot(abs(sv1));
