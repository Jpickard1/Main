
clear; clc;
%% Read in all polymer structures to Kronecker hypergraph
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 200;
EC = cell(nsamps,3);
path2data = "C:\Joshua\MissingData\Projects\KroneckerHypergraphPaper\3CAnalysis\HipHop2HG\";
for s=1:3
    for i=1:nsamps
        filename = "adjList_" + SIMPARMS(s) + "_" + string(i) + "_0.01.txt";
        filepath = path2data + filename; disp(filepath);
        E = readmatrix(filepath);
        EC{i,s} = E;
    end
end

%% Compute clique expansions
CC = cell(size(EC));
for s=1:3
    for i=1:nsamps
        disp(i);
        E = EC{i,s};
        A = zeros(5000, 5000);
        for e=1:size(E)
            hedge = E(e,:);
            A(hedge(1), hedge(2)) = 1;
            A(hedge(1), hedge(3)) = 1;
            A(hedge(2), hedge(1)) = 1;
            A(hedge(2), hedge(3)) = 1;
            A(hedge(3), hedge(1)) = 1;
            A(hedge(3), hedge(2)) = 1;
        end
        CC{i,s} = sparse(A);
    end
end

%% Largest eigenvalue of A
maxLambdaA = zeros(200,3);
for s=1:3
    for i=1:nsamps
        disp(i);
        A = CC{i,s};
        maxLambdaA(i,s) = eigs(A,1);
    end
end

figure;
boxplot(maxLambdaA);
xlabel('Simulation Parameters');
xticklabels(["OFF", "ON", "HIGH"]);
ylabel('Largest Eigenvalue');
title('Eigenvalues of A');

%% Smallest eigenvalue of L
minLambdaL = zeros(200,3);
for s=1:3
    for i=1:nsamps
        disp(i);
        A = CC{i,s};
        D = diag(sum(A));
        L = D - A;
        minLambdaL(i,s) = eigs(L,1,1e-9);
    end
end

figure;
boxplot(minLambdaL);
xlabel('Simulation Parameters');
xticklabels(["OFF", "ON", "HIGH"]);
ylabel('Smallest Eigenvalue');
title('Eigenvalues of L');
