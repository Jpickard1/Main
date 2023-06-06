%% HyperKronFit Driver
%
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023
clear; clc

n = 8;
k = 2;

HG = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) * 0.3), k);
A = HG.adjTensor;

[theta] = HyperKronFit(A)

%% Unit test for edgeProbability

for n0=2:5
    theta = rand(n0,n0);
    A = kron(theta,kron(theta,kron(theta,kron(theta,theta))));
    n = size(A,1);
    disp(n);
    for i=1:n
        for j=1:n
            p1 = edgeProbability(n, theta, i, j);
            assert(A(i,j) == p1);
        end
    end
end
%%
clear; clc;
n0 = 2;
theta = rand(n0,n0); theta = theta ./ sum(theta); theta = theta + theta';
A = kron(theta,kron(theta,kron(theta,kron(theta,theta)))); 
A = (A > median(median(A)));

p=1:32;
j = 11;
k = 1;
p([j k]) = p([k j]);
p2 = permutationProbabilityRatio(p, theta, A, j, k)



%%
B = [1 0;
     0 1];
A = kron(B, B);
p=1:4;
j = 2;
k = 1;
p([j k]) = p([k j]);
permutationProbabilityRatio(p, B, A, j, k)

%% Unit test for samplePermutation
n = 3;
k = 2;
HG = HAT.uniformErdosRenyi(n, round(nchoosek(n,k) * 0.6), k);
A = HG.adjTensor;
C = kron(A, kron(A,A))

p = samplePermutation(C, A)

%% Unit test naive kron fit
A = erdos_renyi_network(64,nchoosek(64,2) * 0.5);
[theta, likelihoods] = NaiveKronFit(A, true, true);
figure; plot(real(likelihoods));

%% Read in snap file
E = readAdjList(filePath, 4);

A = sparse(E(:,1), E(:,2), 1);

[theta, likelihoods] = NaiveKronFit(A, true, true);



%%
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
B = erdos_renyi_network(3,2);
A = kron(B, B);
A = kron(B, A);

%%
B = [1 1;
     1 0];
A = kron(B,B);
p = samplePermutation(A,B)

%%
[theta, l] = NaiveKronFit(real(A), true)
figure; plot(real(l))
