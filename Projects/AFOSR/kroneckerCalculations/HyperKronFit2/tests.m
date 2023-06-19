%% Tests
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 18, 2023

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
[x, y] = find(A == 1); E = [x y];
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