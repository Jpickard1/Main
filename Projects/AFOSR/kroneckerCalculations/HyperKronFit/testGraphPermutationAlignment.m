%% Test Graph Alignment
%
%   Graph alignment is a critical task for implementing the kronfit
%   algorithm. This file serves to test the functions responsible for this
%   implementation.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 7, 2023

%% Test 1
%   Idea: does the algorithm work to correctly align 2 graphs where I know
%   the alignment. We test this with the adjacency matrix that is the 4x4 
%   antidiagonal matrix. The set of node permutations that are ''correct''
%   for such a matrix is all pairs (i, n-i) when there are i vertices for
%   i=1...n. The following test ensures that graphs of this form are
%   aligned correctly.

kronExp = 3;
theta = [0.1 0.9;
         0.9 0.1];
P = theta;
for i=1:kronExp
    P = kron(theta,P);
end
A = (P > 0.5);
p = firstPermutation(A, theta, 1000);
a = length(unique(sum([p; flip(p)],1)));
assert(a == 1);

%% Test 2
%   Idea: watch how the graph becomes aligned over the course of execution

kronExp = 3;
eps = 0.1;
theta = [1 1 1;
         0 1 0
         1 0 1];
theta = theta - eps; theta(theta < 0) = eps;

P = theta;
for i=1:kronExp
    P = kron(theta,P);
end
A = (P > 0.25);
[p, pp] = firstPermutation(A, theta, 50000);

figure; 
subplot(1,3,1); imagesc(A); title('Kronecker Graph');
subplot(1,3,2); imagesc(P); title('Kronecker Expansion');

% Check graph alignment of A into P using the perms
n = size(A,1);
Ap = zeros(n,n);
for i=1:n
    for j=1:n
        Ap(i,j) = A(p(i), p(j));
    end
end
subplot(1,3,3); imagesc(Ap); title('Aligned Graph');

%% Show graph alignment as a movie
h = figure;
h.Visible = 'off';
M(size(pp,1)) = struct('cdata',[],'colormap',[]);
for t = 1:size(pp,1)
    Ap = zeros(n,n);
    for i=1:n
        for j=1:n
            Ap(i,j) = A(pp(t,i), pp(t,j));
        end
    end
    imagesc(Ap);
    drawnow
    M(t) = getframe;
end
% h.Visible = 'on';
figure; movie(M);

