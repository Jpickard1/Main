function [isObsv] = isObsvHG(HG, C)
%ISOBSVHG Checks if a hypergraph and output vertices are observable
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 7, 2023

A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
n = size(A, 1);             % Get n
k = length(size(A));        % Get k

% Unfold A
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

symVars = symvar(sym('x_%d',[n 1]));  % Get symbolic variables 

J = cell(n,1);
P = sparse(Amat);               % Tracks terms AB1B2B... in equation 6
for d=1:n
    if d ~= 1 && d ~= 2
        % getBp is equation 7 in overleaf document 
        P = P * getBp(sparse(Amat), d, k);
        sVec = loadSymVecs(n, size(P,2));
        J{d} = P * sVec;            % Equation 6
    elseif d == 1
        J{1} = sym('x_%d',[n 1]);       % J0 in latex
    elseif d == 2
        sVec = loadSymVecs(n, size(Amat,2));
        J{2} = Amat * sVec;             % J1 in latex
    end
    % Check if the system is observable after each Lie derivative is
    % computed
    Oi = cell(d,1);             % Compute first equality in equation 10
    for i=1:d
        Oi{i} = C * J{i};      % Compute first equality in equation 10
    end
    Oimat = sym([n,n]);         
    for r=1:d                   % Set symbolic matrix to save equation 10 for vertex vx
        for v=1:n               % Compute second equality in equation 10
            Oimat(r,v) = gradient(Oi{r}, symVars(v));
        end
    end
    if rank(Oimat) == n
        isObsv = true;
        return;
    end
end

isObsv = false;

end