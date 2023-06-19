%% Tests
%
% GRAPH ALIGNMENT 1
%   Purpose: Test code's ability to identify a correct mapping between 2
%            graphs using metropolis sampling algorithm
%   Functions:
%    - genPermutation.m
%    - permutationProbabilityRatioTest.m
%
% HYPEREDGE LOG LIKELIHOODS
%   Purpose: A simple verificiation that the hyperedge log likelihoods are
%            computed correctly on small hypergraphs where I can perform
%            these calculations by hand
%   Functions:
%    - hedgeLL.m
%    - noHEdgeLL.m
%
% EMPTY HYPERGRAPH LOG LIKELIHOOD
%   Purpose: verify the log likelihood of empty hypergraphs is computed
%            correctly from theta. This test is for graphs or 2-hypergraphs
%   Functions:
%    - emptyLL.m
%
% KRONECKER HELPER
%   Purpose: verify that the helper function are working as intended
%   Functions:
%    - getCounts.m
%    - kronIndices.m
%
% HYPEREDGE LOG LIKELIHOOD DERIVATIVES
%   Purpose: verify the gradient of the log likelihood is computed
%            correctly
%   Functions:
%    - hedgeDLL.m
%    - noHEdgeDLL.m
%    - emptyDLL.m
%
%    - sampleGradient.m
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

%% HYPEREDGE LOG LIKELIHOOD DERIVATIVES
clear all; close all; clc;
randSeed = 2; rng(randSeed);
n0 = 2;
theta = rand(n0, n0);
kExp = 1; n = n0^kExp;
i = 1; j = 1; DLL = hedgeDLL(n, theta, [i, j]); assert(DLL(i,j) == 1 / theta(i,j));
i = 1; j = 2; DLL = hedgeDLL(n, theta, [i, j]); assert(DLL(i,j) == 1 / theta(i,j));
i = 2; j = 1; DLL = hedgeDLL(n, theta, [i, j]); assert(DLL(i,j) == 1 / theta(i,j));
i = 2; j = 2; DLL = hedgeDLL(n, theta, [i, j]); assert(DLL(i,j) == 1 / theta(i,j));

i = 1; j = 1; DLL = noHEdgeDLL(n, theta, [i, j]); assert(DLL(i,j) == - 1 / (1 - theta(i,j)));
i = 1; j = 2; DLL = noHEdgeDLL(n, theta, [i, j]); assert(DLL(i,j) == - 1 / (1 - theta(i,j)));
i = 2; j = 1; DLL = noHEdgeDLL(n, theta, [i, j]); assert(DLL(i,j) == - 1 / (1 - theta(i,j)));
i = 2; j = 2; DLL = noHEdgeDLL(n, theta, [i, j]); assert(DLL(i,j) == - 1 / (1 - theta(i,j)));

kExp = 8; n = n0^kExp;
DLL = zeros(size(theta));
for i=1:n; for j=1:n
    DLL = DLL + noHEdgeDLL(n, theta, [i j]);
end; end
assert(norm(DLL - emptyDLL(n, theta)) < 1e-2)

% kExp = 3; n = n0^kExp;
% for i=1:n; for j=1:n
%         assert(isequal(hedgeDLL(n, theta, [i, j]), getCounts(n, theta, [i j]) ./ theta))
% end
% end

%% KRONECKER HELPER
clear all; close all; clc;
n0 = 2;
kExp = 3;
n = n0^kExp;
assert(isequal([1 1 1], kronIndices(1, n, n0)));
assert(isequal([1 1 2], kronIndices(2, n, n0)));
assert(isequal([1 2 1], kronIndices(3, n, n0)));
assert(isequal([1 2 2], kronIndices(4, n, n0)));
assert(isequal([2 1 1], kronIndices(5, n, n0)));
assert(isequal([2 1 2], kronIndices(6, n, n0)));
assert(isequal([2 2 1], kronIndices(7, n, n0)));
assert(isequal([2 2 2], kronIndices(8, n, n0)));
n = 4;
assert(isequal([2 0; 0 0], getCounts(n,zeros(2,2),[1,1])));
assert(isequal([1 1; 0 0], getCounts(n,zeros(2,2),[1,2])));
assert(isequal([1 1; 0 0], getCounts(n,zeros(2,2),[1,3])));
assert(isequal([0 2; 0 0], getCounts(n,zeros(2,2),[1,4])));
assert(isequal([1 0; 1 0], getCounts(n,zeros(2,2),[2,1])));
assert(isequal([1 0; 0 1], getCounts(n,zeros(2,2),[2,2])));
assert(isequal([0 1; 1 0], getCounts(n,zeros(2,2),[2,3])));
assert(isequal([0 1; 0 1], getCounts(n,zeros(2,2),[2,4])));
assert(isequal([1 0; 1 0], getCounts(n,zeros(2,2),[3,1])));
assert(isequal([0 1; 1 0], getCounts(n,zeros(2,2),[3,2])));
assert(isequal([1 0; 0 1], getCounts(n,zeros(2,2),[3,3])));
assert(isequal([0 1; 0 1], getCounts(n,zeros(2,2),[3,4])));
assert(isequal([0 0; 2 0], getCounts(n,zeros(2,2),[4,1])));
assert(isequal([0 0; 1 1], getCounts(n,zeros(2,2),[4,2])));
assert(isequal([0 0; 1 1], getCounts(n,zeros(2,2),[4,3])));
assert(isequal([0 0; 0 2], getCounts(n,zeros(2,2),[4,4])));

%% EMPTY HYPERGRAPH LOG LIKELIHOOD
%   An assert tolerance is used because emptyLL is an approximation
clear all; close all; clc;
randSeed = 2; rng(randSeed);
kExp = 5;
n0 = 2;
theta = rand(n0, n0);
ll = 0;
for i=1:n0^kExp; for j=1:n0^kExp
        ll = ll + noHEdgeLL(n0^kExp, theta, [i j]);
end; end
assert(abs((emptyLL(n0^kExp, theta) - ll) / ll) < 1e-3);

kExp = 8;
theta = rand(n0, n0);
ll = 0;
for i=1:n0^kExp; for j=1:n0^kExp
        ll = ll + noHEdgeLL(n0^kExp, theta, [i j]);
end; end
assert(abs((emptyLL(n0^kExp, theta) - ll) / ll) < 1e-3); 

%% HYPEREDGE LOG LIKELIHOODS
clear all; close all; clc;
randSeed = 1; rng(randSeed);
theta = rand(2,2);
n = 2;
assert(hedgeLL(n, theta, [1 1]) == log(theta(1,1)))
assert(hedgeLL(n, theta, [1 2]) == log(theta(1,2)))
assert(hedgeLL(n, theta, [2 1]) == log(theta(2,1)))
assert(hedgeLL(n, theta, [2 2]) == log(theta(2,2)))
assert(noHEdgeLL(n, theta, [1 1]) == log(1 - theta(1,1)))
assert(noHEdgeLL(n, theta, [1 2]) == log(1 - theta(1,2)))
assert(noHEdgeLL(n, theta, [2 1]) == log(1 - theta(2,1)))
assert(noHEdgeLL(n, theta, [2 2]) == log(1 - theta(2,2)))

n = 4;
assert(hedgeLL(n, theta, [1 1]) == log(theta(1,1)) + log(theta(1,1)));
assert(hedgeLL(n, theta, [1 2]) == log(theta(1,1)) + log(theta(1,2)));
assert(hedgeLL(n, theta, [1 3]) == log(theta(1,2)) + log(theta(1,1)));
assert(hedgeLL(n, theta, [1 4]) == log(theta(1,2)) + log(theta(1,2)));
assert(hedgeLL(n, theta, [2 1]) == log(theta(1,1)) + log(theta(2,1)));
assert(hedgeLL(n, theta, [2 2]) == log(theta(1,1)) + log(theta(2,2)));
assert(hedgeLL(n, theta, [2 3]) == log(theta(1,2)) + log(theta(2,1)));
assert(hedgeLL(n, theta, [2 4]) == log(theta(1,2)) + log(theta(2,2)));
assert(hedgeLL(n, theta, [3 1]) == log(theta(2,1)) + log(theta(1,1)));
assert(hedgeLL(n, theta, [3 2]) == log(theta(2,1)) + log(theta(1,2)));
assert(hedgeLL(n, theta, [3 3]) == log(theta(2,2)) + log(theta(1,1)));
assert(hedgeLL(n, theta, [3 4]) == log(theta(2,2)) + log(theta(1,2)));
assert(hedgeLL(n, theta, [4 1]) == log(theta(2,1)) + log(theta(2,1)));
assert(hedgeLL(n, theta, [4 2]) == log(theta(2,1)) + log(theta(2,2)));
assert(hedgeLL(n, theta, [4 3]) == log(theta(2,2)) + log(theta(2,1)));
assert(hedgeLL(n, theta, [4 4]) == log(theta(2,2)) + log(theta(2,2)));
assert(noHEdgeLL(n, theta, [1 1]) == log(1 - exp(log(theta(1,1)) + log(theta(1,1)))));
assert(noHEdgeLL(n, theta, [1 2]) == log(1 - exp(log(theta(1,1)) + log(theta(1,2)))));
assert(noHEdgeLL(n, theta, [1 3]) == log(1 - exp(log(theta(1,2)) + log(theta(1,1)))));
assert(noHEdgeLL(n, theta, [1 4]) == log(1 - exp(log(theta(1,2)) + log(theta(1,2)))));
assert(noHEdgeLL(n, theta, [2 1]) == log(1 - exp(log(theta(1,1)) + log(theta(2,1)))));
assert(noHEdgeLL(n, theta, [2 2]) == log(1 - exp(log(theta(1,1)) + log(theta(2,2)))));
assert(noHEdgeLL(n, theta, [2 3]) == log(1 - exp(log(theta(1,2)) + log(theta(2,1)))));
assert(noHEdgeLL(n, theta, [2 4]) == log(1 - exp(log(theta(1,2)) + log(theta(2,2)))));
assert(noHEdgeLL(n, theta, [3 1]) == log(1 - exp(log(theta(2,1)) + log(theta(1,1)))));
assert(noHEdgeLL(n, theta, [3 2]) == log(1 - exp(log(theta(2,1)) + log(theta(1,2)))));
assert(noHEdgeLL(n, theta, [3 3]) == log(1 - exp(log(theta(2,2)) + log(theta(1,1)))));
assert(noHEdgeLL(n, theta, [3 4]) == log(1 - exp(log(theta(2,2)) + log(theta(1,2)))));
assert(noHEdgeLL(n, theta, [4 1]) == log(1 - exp(log(theta(2,1)) + log(theta(2,1)))));
assert(noHEdgeLL(n, theta, [4 2]) == log(1 - exp(log(theta(2,1)) + log(theta(2,2)))));
assert(noHEdgeLL(n, theta, [4 3]) == log(1 - exp(log(theta(2,2)) + log(theta(2,1)))));
assert(noHEdgeLL(n, theta, [4 4]) == log(1 - exp(log(theta(2,2)) + log(theta(2,2)))));

%% GRAPH ALIGNMENT 1
clear all; close all; clc;
randSeed = 1; rng(randSeed);
kronExp = 3;
eps = 0.1;
theta = [1 1 1;
         0 1 0;
         1 0 1];
theta = theta - eps; theta(theta < 0) = eps;

P = theta;
for i=1:kronExp
    P = kron(theta,P);
end
A = (P > 0.25);
[x, y] = find(A == 1); E = [x y]; n = size(A,1);
tic;
p = randperm(n);
[p, pp] = genPermutation(p, theta, 100000, E);
t1 = toc;

% Check graph alignment of A into P using the perms
n = size(A,1);
Ap = zeros(n,n);
for i=1:n
    for j=1:n
        Ap(i,j) = A(p(i), p(j));
    end
end

figure; 
subplot(1,3,1); imagesc(A); title('Kronecker Graph');
subplot(1,3,2); imagesc(P); title('Kronecker Expansion');
subplot(1,3,3); imagesc(Ap); title('Aligned Graph');

h = figure;
M(size(pp,1)) = struct('cdata',[],'colormap',[]);
for t = 1:size(pp,1)
    Ap = zeros(n,n);
    for i=1:n; for j=1:n
        Ap(i,j) = A(pp(t,i), pp(t,j));
    end; end
    imagesc(Ap);
    M(t) = getframe;
end