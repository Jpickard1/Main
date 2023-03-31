function [Jp] = Jp2Recursive(HG, p, S)
%JP2Recursive
%   
%   HG is a hypergraph from which we extract:
%       Amat: unfolded adjacency tensor
%       k   : order of hypergraph
%       n   : number of vertices in hypergraph
%   p   : which derivative we are interested in
%   S   : sets of vectors of size n
%
%   This function computes the Jp vectors in the shared overleaf document
%   in a reasonable amount of time. At no point is a matrix larger than
%   Amat stored, and there are never more than k-2 kronecker products taken
%   consecutively.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 30, 2023

disp('JP 2 Recursive Called');
disp("    p:" + string(p));

% Extract Amat, n, and k
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix
n = size(Amat,1);
k = length(size(A));

if p == 1
    for i=1:length(S)
        Si = S{i};
        x = Si{1};
        for j=2:length(Si)
            x = kron(x, Si{j});
        end
        if i == 1
            d = Amat * x;
        else
            d = d + Amat * x;
        end
    end
end

% If p ~= 1, then we make the sets recursively and call the function again

d = Jp2Recursive(HG, p-1, Snew);
end

%{
if n < maxP; disp('maxP must be less than n.'); Jp = 0; return; end;
if nargin == 2; x = sym('x_%d',[n 1]); end      % Symbolic state vector
if maxP == 1; Jp = x; return; end;



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
    if i==1
        Jp = Amat * xx;
        % X = xx; 
    else
        X = Amat * xx;
        Jp = Jp + X;
    end
end

% Jp = Amat * X;

end
%}