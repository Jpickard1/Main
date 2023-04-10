function [Jp] = Jp3(HG, p, S)
%JP This function computes the Jp vectors in the shared overleaf document
%   in a reasonable amount of time. At no point is a matrix larger than
%   Amat stored, and there are never more than k-2 kronecker products taken
%   consecutively.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 31, 2023

% Extract Amat, n, k
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix
n = size(Amat,1);
k = length(size(A));

% Recursive base case
if p == k - 2
    xx = S(:,1);
    for i=2:size(S, 2)
        xx = kron(xx, S(:,i));
    end
    Jp = Amat * xx;
    return
end

% Recursive case

% Calculate new sets
b = (p-1)*k-(2*p-3);
disp(b);
Snew = cell(b, 1);
for j=1:length(Snew)
    ss = sym('x_%d', [n, b]);
    offset = 0;
    for i=1:b
        if i ~= j
            ss(:,i) = S(:,i+offset);
        else
            xx = S(:,i);
            offset = k-2;
            for l=1:k-2
                xx = kron(xx, S(:,i+l));
            end
            ss(:,i) = Amat * xx;
        end
    end
    Snew{j} = ss;
end

% Recursively call each set
for i=1:length(Snew)
    if i ~= 1
        Jp = Jp + Jp3(HG, p-1, Snew{i});
    else
        Jp = Jp3(HG, p-1, Snew{i});
    end
end

end

%{
p = maxP;
xInit = p*k - (2*p-1);
b = (p-1)*k-(2*p-3);
S = cell(b, 1);             % S is a cell array of matrices. There are (p-1)*k-(2p-3)
                            % matrices saved in total (i.e. the number of
                            % summations required to construct Bp, and
                            % there are p*k-(2p-1) vectors per matrix,
                            % indicating the individual vectors multiplied
                            % with each matrix term when defining Bp.

for i=1:b
     %ss = zeros(n, b);
    S{i} = repmat(x, 1, xInit);
end

for p=maxP:-1:2   % Loop over Bp Bp-1 ... B3 B2 A
    disp(p)
    b = (p-1)*k-(2*p-3);
    Snew = cell(b,1);
    for j=1:b   % Loop over Si
        if isnumeric(x)
            ss = zeros(n,b);
        else
            ss = sym('x', [n, b]);
        end
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
end
for i=1:size(S,1)
    ss = S{i};
    xx = ss(:,1);
    for j=2:size(ss,2)
        xx = kron(xx, ss(:,j));
    end
    if i==1; X = xx; else X = X + xx; end
end

Jp = Amat * X;

end

%}