% DEMO: Identifying minimal observable nodes in uniform hypergraph
%
%   Overleaf document: https://www.overleaf.com/7146613317rhvrjgqbkprt
%
%   REQUIREMENTS: The following toolboxes are needed to run this script
%       - Symbolic toolbox
%       - Tensor toolbox
%       - Hypergraph Analysis Toolbox
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 17, 2023

%% Hypergraph constructor
clear; close; clc;
% TODO: Set form of k-uniform incidence matrix
N = 5; k = 4;
HG = hyperring(N, k);

% Construct hypergraph
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object

% Unfold A
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

for n=1:N

    % Symbolic vars
    x = sym('x_%d',[N 1]);      % Set symbolic state vector
    symVars = symvar(x);        % Get symbolic variables 
    
    % Compute J vectors (equations 5 and 6 in overleaf document)
    % This is the time computationally expensive part of the program
    J = cell(n,1);
    J{1} = x;                       % J0 in latex
    J{2} = Amat * vecPower(x,k-1);  % J1 in latex
    P = sparse(Amat);               % Tracks terms AB1B2B... in equation 6
    for i=2:n
        % getBp is equation 7 in overleaf document 
        if i~=1
            P = P * getBp(sparse(Amat), i, k);
        else
            P = sparse(Amat);               % Tracks terms AB1B2B... in equation 6
        end
        J{i} = P * vecPower(x, i*k-(2*i-1));  % Equation 6
        % disp(i);                  % Outputs current progress of code
    end
    
    % Compute observability matrices for all vertices
    O = cell(n,1);                  % Cell to hold observability matrices for all vertices
    for vx=1:N
        Ci = zeros(1,N); Ci(vx) = 1;% Equation 9
        Oi = cell(n,1);             % Compute first equality in equation 10
        for i=1:n
            Oi{i} = Ci * J{i};      % Compute first equality in equation 10
        end
        Oimat = sym([n,n]);
        for i=1:n                   % Set symbolic matrix to save equation 10 for vertex vx
            for j=1:N               % Compute second equality in equation 10
%                f = symfun(Oi{i}, [Oi{i}]);
                Oimat(i,j) = gradient(Oi{i}, symVars(j));
            end
        end
        O{vx} = Oimat;              % Save observability matrix for specific vertex
    end
    
    % Greedy Node Selection for Minimum Observable Nodes
    S = 1:N;                                    % Unobserved vertices
    D = [];                                     % Observed vertices
    OD = [];                                    % Observability Matrix
    while rank(OD) < N
        deltaS = zeros(length(S),1);              % Vector to store all changes in rank
        for i=1:length(S)                       % Try all vertices in S
            vx = S(i);                          % Get vertex
            ODS = [OD; O{vx}];                  % Set possible new observabliity matrix
            deltaS(i) = rank(ODS) - rank(OD);      % Compute improved rank for vx
        end
        [~, vx] = max(deltaS);                  % Greedy selection of vertex
        D = [D S(vx)];                          % Add new vertex to observe
        OD = [OD; O{S(vx)}];                    % Set new observability matrix
        S = S(S ~= S(vx));                      % Remove vx from S
    end
    
    disp(OD);   % Display observability matrix
    disp(D);    % Display selected nodes

end