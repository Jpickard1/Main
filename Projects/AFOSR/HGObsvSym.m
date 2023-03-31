function [O] = HGObsvSym(HG)
%HGOBSVSYM This funciton computes symbolic observability matrices for
%   hypergraphs.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 20, 2023
%
% Mod : March 31, 2023 - I updated this file to use Jp3. I then created a
%                        file HGObsvSym0 which used the original Jp.

n = size(HG.IM,1);

% Compute Jp vectors
J = cell(n+1,1);
J{1} = sym('x', [n, 1]);
for i=n:-1:1
    disp(i); p = i; k = 3;
    rpts = p*k-(2*p-1);
    xx = sym('x', [n, 1]);
    S = repmat(xx, 1, rpts);
    J{i+1} = Jp3(HG, i, S);
end

% Compute observability matrices for all vertices
O = cell(n,1);                  % Cell to hold observability matrices for all vertices
symVars = symvar(sym('x',[n 1]));  % Get symbolic variables 
for vx=1:n
    disp(vx);
    Ci = zeros(1,n); Ci(vx) = 1;% Equation 9
    Oi = cell(n+1,1);             % Compute first equality in equation 10
    for i=1:n+1
        Oi{i} = Ci * J{i};      % Compute first equality in equation 10
    end
    Oimat = sym([n+1,n]);         % + 1 added since we need 0,...,n
    for i=1:n+1                   % Set symbolic matrix to save equation 10 for vertex vx
        for j=1:n               % Compute second equality in equation 10
            Oimat(i,j) = gradient(Oi{i}, symVars(j));
        end
    end
    O{vx} = Oimat;              % Save observability matrix for specific vertex
end

end
