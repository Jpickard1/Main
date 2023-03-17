% Identifying observable nodes in uniform hypergraph
%
%   DOC: https://www.overleaf.com/7146613317rhvrjgqbkprt
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 12, 2023

%% 03/17/2023 - Tensor Eigenvectors
clear

T(:,:,1) = [1 2; 2 3]; T(:,:,2)= [2 3; 3 6];
T = tensor(T)
issymmetric(T)

IM = [1 1 1;
      1 1 0;
      1 0 1;
      0 1 1];
k = 3;
[n,e] = size(IM);
HG = Hypergraph('IM', IM);  % Hypergraph object
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);   

[l, v] = heig(A)

A = tensor(A);   
A11 = ttv(A, v(:,1), 1)
A12 = ttv(A11, v(:,1), 1)

A21 = ttv(A, v(:,2), 1)
A22 = ttv(A21, v(:,2), 1)


%% 03/06/2023 - Are LTI Systems observable up to a scaler output?

clear

A = rand(2,2);
C = [1 0];

lqe(A,zeros(size(A)),C,zeros(size(A)),ones([1 1]),zeros([2 1]))

%% 03/05/2023 - Random calculations to verify kron ttv kron == ttv kron ttv
i = 0;
F = [];
while isempty(F)
    n1 = randi([1,10],1);
    n2 = randi([1,10],1);
    A = rand(n1,n1,n1);
    B = rand(n2,n2,n2);
    C = rand(n1,1);
    D = rand(n2,1);
    
    AC = tenmat(ttv(tensor(A),C,1),1); AC = AC(:,:);
    BD = tenmat(ttv(tensor(B),D,1),1); BD = BD(:,:);
    ACBD = kron(AC,BD);
    
    AB = superkron(A,B);
    CD = kron(C,D);
    ABCD = tenmat(ttv(tensor(AB), CD, 1), 1); ABCD= ABCD(:,:);
    
    % Check number of different values
    F = find(((ABCD - ACBD) < 1e-10) == 0);
    i = i + 1; disp(i);
end

% Final Output
%
% >> 3286156

%% 03/05/2023 - Symbolic calculations at Elises
clear; close all; clc

A = sym('a', [2,2,2]);
B = sym('b', [2,2,2]);
x = sym('x', [2, 1]);
y = sym('y', [2, 1]);

Ax = ttv(tensor(A),x,1);
By = ttv(tensor(B),y,1);

Axm = sym('x', size(Ax)); s1 = size(Ax); n1=s1(1); m1=s1(2);
for i=1:n1; for j=1:m1; Axm(i,j) = Ax(i,j); end; end
Bym = sym('x', size(By)); s2 = size(By); n2=s2(1); m2=s2(2);
for i=1:n2; for j=1:m2; Bym(i,j) = By(i,j); end; end

Om = kron(Axm,Bym)

%% 03/04/2023
clear; close all; clc

V = 10; k=3;
HG = hyperring(V, k);
A = HG.adjTensor;
A = tensor(A);              % A tensor as tensor object

% Unfold A
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

p=4;

% START: getp.m function
n = size(Amat, 1);
B = sparse(n^(p+1), size(Amat, 2) * n^(p+1));    % Set sparse matrix to return
bound = (p-1)*k-(2*p-3);        % Compute upper bounds of loop
for i=1:bound
    if i==1
        P = sparse(Amat);
    else
        P = sparse(eye(n,n));
    end
    for j=2:bound
        if j ~= i
            P = kron(P, eye(n,n));
        else
            P = kron(P, Amat);
        end
    end
    if i~=1
        B = B + P;
    else
        B = P;
    end
    disp("=====");
    disp(rank(full(P)));
    disp(rank(full(B)));
end
% END: getp.m function

%% Kronecker sum for rank

clear;

V = 6;
A = rand(V,2 * V);
A(V,:) = A(1,:);

Ak = kronSum(A,A);
disp(rank(A));
disp(rank(Ak));

Iv = eye(V);
S1 = kron(kron(kron(Iv,A), Iv), Iv); % rank = rA * n^3 = 625
S2 = kron(kron(kron(Iv,Iv), A), Iv); % rank = rA * n^3 = 625
S = S1 + S2;

B = S1/S2;

rank(S1)
rank(S2)
rank(S)

% B = rand(V,V);
% I = eye(V,V);
% kronSum = kron(A,I) + kron(I,B);
% rank(kronSum)


%% 03/02/2023 - night
clear; close all; clc

n1 = 2;
n2 = 3;
p = 1;
q = 2;

A1 = sym('a_%d%d', [n1, n1]);
A2 = sym('a_%d%d', [n2, n2]);
C1 = sym('c_%d%d', [p, n1]);
C2 = sym('c_%d%d', [q, n2]);
Ak = kron(A1, A2);
Ck = kron(C1, C2);

O1 = obsvSym(A1, C1);
O2 = obsvSym(A2, C2);
Ok = obsvSym(Ak, Ck)


%% 03/02/2023
clear; close all; clc

n1 = 5;
n2 = 6;

A1 = eye(n1);
A2 = eye(n2);
C1 = rand(2,n1);
C2 = rand(2,n2);
Ak = kron(A1, A2);
Ck = kron(C1, C2);

O1 = obsv(A1, C1);
O2 = obsv(A2, C2);
Ok = obsv(Ak,Ck);

O12 = kron(O1,O2);

CC = (Ok == O12);

% Compare O12 and Ok
for i=1:size(Ok,1)
    rowi = Ok(i,:);
    for j=1:size(O12,1)
        rowj = O12(j,:);
        if sum(rowi == rowj) == size(Ok,2)
            disp(i)
        end
    end
end





%% 2/17-19/2023
D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');
[M, idxs] = HAT.multicorrelations(D, 3, 'Wang');
figure; histogram(M)

t = 0.83; % Disconnected at 85
hyperedges = idxs(find(M>t),:);
IM = HAT.hyperedge2IM(hyperedges);
HG = Hypergraph('IM', IM);
figure; plot(graph(HG.cliqueGraph)) % To check the graph is fully connected

max(M)
sum(M>0.85)

%% Hypergraph constructor
%{
IM = [1 1 1;                % Hyperstar incidence matrix
      1 0 0;
      1 0 0;
      0 1 0;
      0 1 0;
      0 0 1;
      0 0 1]
%}
IM = [1 1 1;
      1 1 0;
      1 0 1;
      0 1 1];
k = 3;
[n,e] = size(IM);
HG = Hypergraph('IM', IM);  % Hypergraph object
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
%% 02/14/2023 @ 8 PM
Amat = tenmat(A,1); % Unfold A
Amat = Amat(:,:);

x = sym('x_%d',[n 1]);         % Symbolic state vector

J = cell(n,1);
J{1} = x;
J{2} = Amat * vecPower(x,1);
P = sparse(Amat);
for i=2:n
    P = P * getBp(Amat, i, k);
    J{i} = P * vecPower(x, i);
end

Jfun = cell(n,1);
for i=1:n
    Jfun{i} = symfun(J{i}, x);
end

%%
jj = J(2)
gradient(jj{1,1}, symvar(jj))

%% Gradients
% fun = Jfun{1}
gradient(fun)
gradient(Jfun{i}, symvar(Jfun{i}))
gradient(J{i}, symvar(J{i}))

gradient(f,symvar(f))
[a b c d] = gradient(ff,symvar(ff))

ff{1}

%%
x = cell(3, 10);
for i = 1:10
    for j = 1:3
        x{j,i} = sprintf('x%d%d',i,j);
    end
end
x = x(:); % Now x is a 30-by-1 vector
x = sym(x, 'real');

B = [1,1,1;-1,1,1;1,-1,1;-1,-1,1];

A = zeros(40,30);
for i=1:10
    A(4*i-3:4*i,3*i-2:3*i) = B;
end

b = zeros(40,1);

c = sym(zeros(1,10));
i = 1:10;
c = (x(3*i-2).^2 + x(3*i-1).^2 + (x(3*i)+1).^2 - 1).';

gradient(c(1), symvar(c(1)))


%%
syms x y z

f(x,y,z) = [2*y*z*sin(x) + 3*x*sin(z)*cos(y); 2*y*z*sin(x) + 3*x*sin(z)*cos(y)]

gradient(f,symvar(f))


%% 2/14/2023 6 PM
Amat = tenmat(A,1); % Unfold A
Amat = Amat(:,:);

x = sym('x_%d',[n 1]);         % Symbolic state vector
J0 = x;                     % J0
J1 = Amat * vecPower(x,1);
J2 = Amat * getBp(Amat, 2, 3) * vecPower(x,2);
J3 = Amat * getBp(Amat, 3, 3) * vecPower(x,3);

%%
J = cell(n,1);
J{1} = x;
J{2} = Amat * vecPower(x,1);
P = Amat;
for i=2:n
    P = P * getBp(Amat, i, k);
    J{i} = P * vecPower(x, i);
end


%% Symbolic computations
%x = [];
%ii = ['x1','x2','x3','x4','x5','x6','x7'];
%for i=1:7
%    x(i) = sym(ii(i))
%end
x = sym('x_%d',[n 1]);         % Symbolic state vector

% syms x [7 1] matrix

J0 = x;                     % J0
J1 = ttvk(A,x);             % J1
% J2 = ttv(J1,x,1);           % 
% J1t = J1(:)

J0fun = symfun(J0, x)
J1fun = symfun(J1, x)

% gradient(J0func, x)
% for i=1:20;  eval(['syms x' num2str(i) ]); end


