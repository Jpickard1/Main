%% Kronecker Calculations
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 31, 2023

%% April 22, 2023 Laplacians

clear all; clc; close all;
n = 4; k = 3;
HG1 = getToyHG(n, k, 'complete'); L1 = HG1.laplacianTensor; A1 = HG1.adjTensor;
HG2 = getToyHG(n, k, 'complete'); L2 = HG2.laplacianTensor; A2 = HG2.adjTensor;
A = superkron(A1, A2); IM = HAT.A32IM(A); HG = Hypergraph('IM',IM);
L = HG.laplacianTensor;
L12 = superkron(L1, L2);

%%
clear; clc; close all;
r = 3; n = 2;
Ac = cell(r,1);
Ac{1} = rand(n,r);
Ac{2} = rand(n,r);
Ac{3} = rand(n,r);
Bc = cell(r,1);
Bc{1} = rand(n,r);
Bc{2} = rand(n,r);
Bc{3} = rand(n,r);

A = zeros(n,n,n); B = zeros(n,n,n);
for i=1:r
    R = outerProduct(reshape(outerProduct(Ac{1}(:,i),Ac{2}(:,i)), [n,n]),Ac{3}(:,i));
    A = A + R;
end
for i=1:r
    R = outerProduct(reshape(outerProduct(Bc{1}(:,i),Bc{2}(:,i)), [n,n]),Bc{3}(:,i));
    B = B + R;
end

AB = superkron(A,B);
Cc = cell(r,1);
for i=1:r
    Cc{i} = kron(Ac{i}, Bc{i});
end
C = zeros(size(AB));
for i=1:r^2
    R = outerProduct(reshape(outerProduct(Cc{1}(:,i),Cc{2}(:,i)), [n^2,n^2]),Cc{3}(:,i));
    C = C + R;
end

E = (C - AB);
sum(E(:))






%%
clear

U1 = [1 3; 2 4];
U2 = [5 7; 6 8];
U3 = [9 11; 10 12];

O1 = outerProduct(outerProduct(U1, U2), U3);
O2 = outerProduct(outerProduct(U1(:,1), U2(:,1)), U3(:,1)) + outerProduct(outerProduct(U1(:,2), U2(:,2)), U3(:,2));

M = max(O2);
while numel(M) ~= 1
    M = max(M);
end
disp(M)

max(O2)

% ttv(tensor(U1(:,1)), U2(:,1), 2)

%% April 11, 2023 Kronecker producting fun hypergraphs

n = 7; k = 4;
HG1 = getToyHG(n, k, 'complete');
HG2 = getToyHG(n, k, 'complete');

A = superkron(HG1.adjTensor, HG2.adjTensor);
HG = Hypergraph('IM', HAT.A42IM(A));
figure; HG.plot();

size(HG.IM)

%% April 11, 2023 Calculating Kronecker Factorization of Toy Hypergraphs

clear; clc; close all;

n = 49;
k=3;

figure;
HG = hyperchain(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,1); HG.plot(); title('Hyperchain');
        subplot(3,3,2); HG1.plot(); title('Factor 1');
        subplot(3,3,3); HG2.plot(); title('Factor 2');
HG = hyperring(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,4); HG.plot(); title('Hyperring');
        subplot(3,3,5); HG1.plot(); title('Factor 1');
        subplot(3,3,6); HG2.plot(); title('Factor 2');
HG = hyperstar(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,7); HG.plot(); title('Hyperstar');
        subplot(3,3,8); HG1.plot(); title('Factor 1');
        subplot(3,3,9); HG2.plot(); title('Factor 2');

saveas(gcf, "Kronecker Factors of Toy Hypergraphs.png");

%% Adjacency tensor to IM

idxs = find(B ~= 0);
[x,y,z] = ind2sub([10, 10, 10],idxs);
E = [x y z]
IM = HAT.hyperedges2IM(E);
HG = Hypergraph('IM', IM);
HG.plot()


HG = Hypergraph('IM', IM);
HG.plot()

%% April 11, 2023 Calculating the degree sequence of Kronecker Graphs

HG = hyperchain(100,3);
A = HG.adjTensor; A = tensor(A);
v = ones(100,1);
a1 = ttv(A,v, 1);
a2 = ttv(a1,v, 1);
% HG.plot()
s1 = sum(full(HG.IM), 2);
m1 = double(a2)
sum(s1 ~= m1)

%% Hamming Hypergraph Similarity Measure
clear; clc;
a = rand(3,3,3);
b = rand(3,3,3);
c = rand(3,3,3);
d = rand(3,3,3);

sum(sum(sum(a)))*sum(sum(sum(b)))
sum(sum(sum(superkron(a,b))))


sum(sum(sum(a)));
sum(sum(sum(d)));
sum(sum(sum(d)));

%% Reshaping sum of 2 tensors vs sum of 2 reshaped tensors
X = rand(3,3,3);
D = rand(3,3,3);
XD = X + D;

xdk = reshape(XD, 3, 9)
x = reshape(X, 3, 9); d = reshape(D, 3, 9);
xd = x + d;

xdk == xd

%% Entropy of kronecker graph
% It seems that the entropy of a hypergraph is not bound by either the sum
% or the product of the entropy of its kronecker factors
clear; clc;
itrs = 10000;
for i=1:itrs
    n1 = randi(100,1); n2 = randi(100,1);
    x = rand(n1, 1); x = x ./ sum(x);
    y = rand(n2, 1); y = y ./ sum(y);
    xy = kron(x,y);
    ex = entropy(x);
    ey = entropy(y);
    exy = entropy(xy);
    if ex + ey < exy
        disp('FAIL');
    end
end

%% Krnocker Factorization
clear; clc;
n = 3;
x = rand(n,1); y = rand(n,1);

z = kronSum(x, y, 4)

% x = sym("x", [n, 1]); y = sym("y", [n, 1]);
% z = kron(x,y) + kron(y, x);

Z = reshape(z, 3, 3);

[a,b,c] = svd(Z)

b(1,1) * kron(a(1,:),a(1,:)) + b(2,2) * kron(a(2,:),a(2,:)) + b(3,3) * kron(a(3,:),a(3,:))
z'

Z - reshape(b(1,1) * (a(:,1) * c(:,1)'), 3, 3)
a
c
Z


b(1,1) * (a(1,:)' * a(1,:)) + b(2,2) * (a(2,:)' * a(2,:)) + b(3,3) * (a(3,:)' * a(3,:))
Z - (b(1,1) * (a(:,1) * a(:,1)')) + (b(2,2) * (a(:,2) * a(:,2)')) + (b(3,3) * (a(:,3) * a(:,3)'))
Z

a * b * c'

