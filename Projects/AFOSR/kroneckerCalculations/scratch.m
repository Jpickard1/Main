%% Kronecker Calculations
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 31, 2023

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

