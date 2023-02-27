%% Symbolic Computations
%
%   The purpose of this file is to precompute the symbolic vectors required
%   to select the minimal observable node sets. Constructing large symbolic
%   vectors is expensive, but their multiplication with sparse matrices is
%   relatively cheap. To improve the speed of analyzing many small
%   hypergraphs (i.e. n<=7) I precompute all possible required symbolic
%   vectors in this script.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 22, 2023

clc; clear; close all;
maxN = 7;
K = 3:5;

for n=3:3
    % Set symbolic variables for HG with n vxc
    x = sym('x_%d',[n 1]);      % Set symbolic state vector
    symVec = sym('x_%d',[n 1]);      % Set symbolic state vector
    % symVars = symvar(x);        % Get symbolic variables 
    % for ki=1:length(K)
    %     k = K(ki);
    %     if k > n
    %         continue;
    %     end
        i = 0;
        while length(symVec) < 1000
        % for i=0:n
            if i~= 0
                disp(i);
                symVec = kronSymVec(symVec, x, i-1);
            end
            disp(length(symVec));
            fileName = "symVecs/" + string(n) + "_" + string(length(symVec)) + ".mat";
            cmd = "save " + fileName + " " + "symVec";
            eval(cmd);
            disp(i); i = i + 1;
        end 
    % end 
end


%% Pre 02/26/2023
% symVectors = cell(maxN, 1);
symVectors = containers.Map;

for n=4:maxN
    % Set symbolic variables for HG with n vxc
    x = sym('x_%d',[n 1]);      % Set symbolic state vector
    symVars = symvar(x);        % Get symbolic variables 
    
    symVectorsN = cell(length(K), 1);
    for ki=1:length(K)
        k = K(ki);

        if k > n
            continue;
        end

        symVectorsNK = cell(length(n), 1);
        symVectorsNK{1} = x;
        disp(length(x));
        symVec = kronSymVec(x, x, k-1);
        symVectorsN{2} = symVec;
        disp(length(symVec));
        for i=3:n
            symVec = kronSymVec(symVec, x, k-1);
            % symVectorsN{i-1} = symVec;
            disp(length(symVec));
        end 
        
        symVectorsN{ki} = symVectorsNK;
    end 
    symVectors(string(n)) = symVectorsN;
end 

%% kronSymVec

function symVec = kronSymVec(symVec, x, pow)
    n = length(x);
    totalSize = n^pow;
    while length(symVec) < totalSize
        symVec = kron(symVec, x);
    end
end