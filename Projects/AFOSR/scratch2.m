clear; close; clc;

%% 03/20/2023 - Working with JpSym
n = 6;
k = 3;

HG = hyperstar(n,k); % disp(full(HG.IM))
% O = HGObsvSym(HG)
x = rand(n,1);
O = HGObsvNum(HG, x)

%% HGObsvSym Function
% Compute Jp vectors
Jp = cell(n,1);
for i=n:-1:1
    disp(i)
    Jp{i} = JpSym(HG, i); %, rand(n,1))
end

% Compute observability matrices for all vertices
O = cell(n,1);                  % Cell to hold observability matrices for all vertices
symVars = symvar(sym('x_%d',[n 1]));  % Get symbolic variables 
for vx=1:n
    Ci = zeros(1,n); Ci(vx) = 1;% Equation 9
    Oi = cell(n,1);             % Compute first equality in equation 10
    for i=1:n
        Oi{i} = Ci * Jp{i};      % Compute first equality in equation 10
    end
    Oimat = sym([n,n]);         
    for i=1:n                   % Set symbolic matrix to save equation 10 for vertex vx
        for j=1:n               % Compute second equality in equation 10
            Oimat(i,j) = gradient(Oi{i}, symVars(j));
        end
    end
    O{vx} = Oimat;              % Save observability matrix for specific vertex
end


%% Cast symbolic expressions to symbolic functions
symVars = symvar(sym('x_%d',[1 n]));
symE = Jp{1}
f(symVars) = symE
gradient(f,symVars)

gradient(symE, symVars)

%% 03/19/2023 - Recursive strategy to computing derivatives of hypergraph dynamics
n = 5;
k = 5;

HG = hyperring(n,k); % disp(full(HG.IM))

A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

maxP = n;
p = maxP;
xInit = p*k - (2*p-1);
b = (p-1)*k-(2*p-3);
S = cell(b, 1);
x = rand(n,1);
% x = sym('x_%d',[n 1]);         % Symbolic state vector
for i=1:b
    ss = zeros(n, b);
    S{i} = repmat(x, 1, xInit);
end
disp('==============');
disp(size(S));
disp(size(S{1}));

for p=maxP:-1:2   % Loop over Bp Bp-1 ... B3 B2 A
    b = (p-1)*k-(2*p-3);
    Snew = cell(b,1);
    for j=1:b   % Loop over Si
        ss = sym('x', [n, b]);
        offset = 0;
        for i=1:b
            if i ~= j
                ss(:,i) = S{j}(:,i+offset);
            else
                xx = S{j}(:,i);
                offset = k-2;
                for l=1:k-2
                    xx = kron(xx, S{j}(:,i+l));
                end
                ss(:,i) = Amat * xx;
            end
            % disp("            " + string(j));
        end
        Snew{j} = ss;
        % disp("        " + string(j));
    end
    S = Snew;
    disp('==============');
    disp(p);
    disp(b);
    disp(size(S));
    disp(size(S{1}));
    % disp("    " + string(p));
end
% disp(size(S)); disp(size(S{1}));
for i=1:size(S,1)
    ss = S{i};
    xx = ss(:,1);
    for j=2:size(ss,2)
        xx = kron(xx, ss(:,j));
    end
    if i==1; X = xx; else X = X + xx; end
end
Jp = Amat * X;
%% Redone above
for p=maxP:-1:2   % Loop over Bp Bp-1 ... B3 B2 A
    b = (p-1)*k-(2*p-1);
    Snew = cell(b,1);
    for i=1:b   % Loop over Si
        ss = sym('x', [n, b]);
        offset = 0;
        for j=1:b
            if j ~= i
                ss(:,j) = S{i}(:,j+offset);
            else
                xx = S{i}(:,j);
                offset = k-2;
                for l=1:k-2
                    xx = kron(xx, S{i}(:,j+l));
                end
                ss(:,j) = Amat * xx;
            end
            disp("            " + string(j));
        end
        Snew{i} = ss;
        disp("        " + string(i));
    end
    S = Snew;
    disp("    " + string(p));
end
disp(size(S)); disp(size(S{1}));