%% HyperKronFit Driver
%
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

n = 8;
k = 2;

HG = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) * 0.3), k);
A = HG.adjTensor;

[theta] = HyperKronFit(A)

%% Unit test for edgeProbability

B = [0.5  0.25;
     0.75 0.2];
A = kron(B,B);
for i=1:4
    for j=1:4
        p1 = edgeProbability(4, B, i, j);
        if A(i,j) ~= p1
            disp('FAIL');
            disp(A(i,j));
            disp(p1);
            disp(i); disp(j)
        end
    end
end

%% Unit test for samplePermutation
k = 5;
B = rand(2,2);
A = kron(B, B);
for i=1:k
    A = kron(B, A);
end

c = 0;
for i=1:1000
    p = samplePermutation(A, B);
    if p(1) == 1; c = c + 1; end
end

%% Evaluate gradient
B = erdos_renyi_network(8, 8);
A = kron(B, B);
A = kron(B, A);
%%
theta = NaiveKronFit(A)

