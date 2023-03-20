clear; close; clc;
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