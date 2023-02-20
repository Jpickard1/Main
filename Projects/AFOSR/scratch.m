% Identifying observable nodes in uniform hypergraph
%
%   DOC: https://www.overleaf.com/7146613317rhvrjgqbkprt
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 12, 2023

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


