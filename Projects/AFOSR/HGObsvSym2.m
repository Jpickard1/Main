function [O] = HGObsvSym2(HG)
%HGOBSVSYM This funciton computes symbolic observability matrices for
%   hypergraphs.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 20, 2023

n = size(HG.IM,1);

% Compute Jp vectors
J = cell(n,1);
for i=n:-1:1
    disp(i)
    J{i} = Jp(HG, i);
end

% Compute observability matrices for all vertices
O = cell(n,1);                  % Cell to hold observability matrices for all vertices
symVars = symvar(sym('x_%d',[n 1]));  % Get symbolic variables 
for vx=1:n
    Ci = zeros(1,n); Ci(vx) = 1;% Equation 9
    Oi = cell(n,1);             % Compute first equality in equation 10
    for i=1:n
        Oi{i} = Ci * J{i};      % Compute first equality in equation 10
    end
    Oimat = sym([n,n]);         
    for i=1:n                   % Set symbolic matrix to save equation 10 for vertex vx
        for j=1:n               % Compute second equality in equation 10
            Oimat(i,j) = gradient(Oi{i}, symVars(j));
        end
    end
    O{vx} = Oimat;              % Save observability matrix for specific vertex
end

end