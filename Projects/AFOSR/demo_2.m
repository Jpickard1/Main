% DEMO: Identifying minimal observable nodes in uniform hypergraph. This
% demo uses precomputed symbolic vectors to improve the speed of
% calculating the MON set.
%
%   Overleaf document: https://www.overleaf.com/7146613317rhvrjgqbkprt
%
%   REQUIREMENTS: The following toolboxes are needed to run this script
%       - Symbolic toolbox
%       - Tensor toolbox
%       - Hypergraph Analysis Toolbox
%       - precomputed symbolic vectors according to Joshua's specific file
%       structure
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 26, 2023

%% Hypergraph constructor
clear; close; clc;
% TODO: Set form of k-uniform incidence matrix
n = 4; k = 3;
HG = hyperring(n, k);

% Construct hypergraph
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object

% Unfold A
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

% Symbolic vars
x = sym('x_%d',[n 1]);      % Set symbolic state vector
symVars = symvar(x);        % Get symbolic variables 

% Compute J vectors (equations 5 and 6 in overleaf document)
% This is the time computationally expensive part of the program
J = cell(n,1);
J{1} = x;                       % J0 in latex
sVec = loadSymVecs(n, size(Amat,2));
J{2} = Amat * sVec;             % J1 in latex
P = sparse(Amat);               % Tracks terms AB1B2B... in equation 6
for i=2:n
    % getBp is equation 7 in overleaf document 
    P = P * getBp(sparse(Amat), i, k);
    sVec = loadSymVecs(n, size(P,2));
    J{i} = P * sVec;            % Equation 6
end

% Compute observability matrices for all vertices
O = cell(n,1);                  % Cell to hold observability matrices for all vertices
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

% Greedy Node Selection for Minimum Observable Nodes
S = 1:n;                                    % Unobserved vertices
D = [];                                     % Observed vertices
OD = [];                                    % Observability Matrix
while rank(OD) < n
    deltaS = zeros(length(S));              % Vector to store all changes in rank
    for i=1:length(S)                       % Try all vertices in S
        vx = S(i);                          % Get vertex
        ODS = [OD; O{vx}];                  % Set possible new observabliity matrix
        deltaS = rank(ODS) - rank(OD);      % Compute improved rank for vx
    end
    [~, vx] = max(deltaS);                  % Greedy selection of vertex
    D = [D S(vx)];                          % Add new vertex to observe
    OD = [OD; O{S(vx)}];                    % Set new observability matrix
    S = S(S ~= vx);                         % Remove vx from S
end

disp(OD);   % Display observability matrix
disp(D);    % Display selected nodes

function symVec = loadSymVecs(n, l)
%LOADSYMVECS This function loads symbolic vectors from .mat files.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 26, 2023
    fName = "symVecs/" + string(n) + "_" + string(l) + string(".mat");
    load(fName)
end

%% 
%{
IM = [1 1;
      1 0;
      1 0;
      0 1;
      0 1];
IM = [1 1 1;
      1 1 0;
      1 0 1;
      0 1 1];
%}

