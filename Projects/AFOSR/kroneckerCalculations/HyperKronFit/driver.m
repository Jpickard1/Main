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
A = erdos_renyi_network(64,round(nchoosek(64,2) * 0.1));
[theta, likelihoods] = NaiveKronFit(A, true, true);
figure; imagesc(A)
figure; plot(real(likelihoods));

%% Unit test naive kron fit on a kronecker hypergraph
kronExp = 4;
theta = [0.9 0.7
         0.6 0.8]; theta = theta / sum(theta, 'all');
E = kronGen(theta, kronExp, 4 * size(theta,1)^kronExp);
A = sparse(size(theta,1)^(kronExp+1), size(theta,1)^(kronExp+1));
A(E(:,1), E(:,2)) = 1;
figure; imagesc(A)


%% Read in snap file
filePath = "C:\Joshua/Software/snap/examples/as20graph.txt"
E = readAdjList(filePath, 4);
A = sparse(E(:,1), E(:,2), 1);
[theta, likelihoods] = NaiveKronFit(A, true, true);

%% kronecker expansion

kronExp = 4;
eps = 0.1;
theta = [1 1;
         0 1];
theta = rand(3,3);
% theta = theta - eps; theta(theta < 0) = eps;

P = theta;
for i=1:kronExp
    P = kron(theta,P);
end
A = (P > 0.03);
% [p, pp] = firstPermutation(A, theta, 50000);

[theta, likelihoods] = NaiveKronFit(A, true, true, 3);

LP = theta;
for i=1:kronExp
    LP = kron(theta,LP);
end

figure; 
subplot(1,3,1); imagesc(A); title('Kronecker Graph');
subplot(1,3,2); imagesc(P); title('Kronecker Expansion');
subplot(1,3,3); imagesc(LP); title('Learned Kronecker Expansion');


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

%% Code to fix kronGen

clear; clc; close all;
theta = rand(2,2);
n0 = size(theta,1);
theta = theta / sum(sum(theta));
kronExp = 10000;
disp(reshape(theta, [1, numel(theta)]))

falls = randsrc(1,kronExp,[1:numel(theta); reshape(theta, [1, numel(theta)])]);
fallIJ = zeros(numel(falls), 2);
for f=1:numel(falls)
    [fallIJ(f,1), fallIJ(f,2)] = ind2sub(size(theta), falls(f));
end

C = zeros(size(theta));
for f=1:numel(falls)
    C(fallIJ(f,1), fallIJ(f,2)) = C(fallIJ(f,1), fallIJ(f,2)) + 1;
end
disp(C / kronExp);
disp(theta)
























%% Debugging KronFit
% 
%   Known Problems:
%       1. log likelihood evaluates as a positive quantity
%       2. algorithm converges to the wrong values
%
% Auth: Joshua Pickard
% Date: June 9, 2023

%% 1. log likelihood evaluates as a positive quantity
%
%   Possible causes - I don't enforce the values of theta to be a
%   probability distribution i.e. sum(sum(theta)) ~= 0 necessarily.
%       Potential fix - I normalize theta prior to evaluating the log
%       likelihood and gradients.
%           This does fix the `symptom` but it is possible that there is 
%           a deeper cause of the issue.


[thetaLearned3, likelihoods3] = NaiveKronFit(A, true, true, 3, thetaLearned2);

%% 2. algorithm converges to the wrong values
%
%   Possible causes: it doesn't evaluate the likelihood or gradient
%   correctly.
%       Potential fix - I am going to normalize theta before computing 
%       the gradient of the log likelihood
%
%   Strategies and thorughts for checking issues in the convergence:
%       - It does appear the algorithm converges, just to the wrong 
%         value of theta
%       - I can try inputting the generator theta and checking how far
%         the code modifies it
%       - I can try putting in theta far from the generator and checking
%         for how it causes the algorithm to converge
%       - I can try rotations/flips of theta
%       - I can try theta near the generating theta
%       - I can try theta0 = ones or theta0 = zeros
%       - I could work an example by hand
%       

clear; clc;
% Make a graph
theta = [1 1 1;         1 1 0;         1 0 1];
% theta = [1 1 0; 1 1 0; 0 0 0];
eps = 1e-3; theta = theta - eps; theta(theta < 0) = eps;
n0 = size(theta, 1); kronExp = 5; numE = 50 * n0^kronExp;
theta = theta / sum(sum(theta));
E = kronGen(theta, kronExp, numE); n = n0^kronExp;
A = sparse(n, n);
for i=1:size(E,1); A(E(i,1), E(i,2)) = 1; end
A = full(A);
figure; imagesc(A)

% theta0 = ones(3,3);
theta0 = rand(3,3);
[thetaLearned, likelihoods, thetas] = NaiveKronFit(A, true, true, 3, theta0, 100);


thetaLearned / sum(sum(thetaLearned))
theta


%% Example by hand
clear; close all; clc
eps = 1e-2;
thetaG = [1-eps 1-eps; 1-eps eps];
A = kron(thetaG,kron(thetaG,kron(thetaG, kron(thetaG, thetaG))));
A = (A>0.5);
figure; imagesc(A);
theta0 = 0.5 * ones(2,2);
theta0 = rand(2,2);
[thetaLearned, likelihoods, thetas] = NaiveKronFit(A, true, true, 3, theta0, 300);
% [thetaLearned2, likelihoods, thetas] = NaiveKronFit(A, true, true, 3, thetaLearned, 100);
% theta0 = sym('t_%d%d',[2 2]);

%% Check getEmptyGrapGrad and getNoEdgeDLL2 are equivalent
theta = theta0;

R1 = getEmptyGraphGrad(n,theta);
R2 = zeros(size(R1));
for u=1:n
    for v=1:n
        theta = theta / sum(sum(theta));
        n0 = size(theta,1);
        kronExp = log(n) / log(n0);
        % edgeP = edgeProbability(n, theta, u, v);
        eLL = edgeLL(n, theta, u, v);
        noEdgeLL = log(1 - exp(eLL));
        % Count the number of times an entry of theta is used
        count = zeros(size(theta));
        for i=1:kronExp
            i1 = mod(floor((u - 1)/ n0^(i-1)), n0) + 1;
            i2 = mod(floor((v - 1)/ n0^(i-1)), n0) + 1;
            count(i1, i2) = count(i1, i2) + 1;
        end
        gradient = zeros(size(theta));
        % Calculate gradient of log likelihood function
        for i=1:size(theta,1)
            for j=1:size(theta,2)
                % Count the number of times (i,j) was used
                c = count(i,j);
                negGrad = getNoEdgeDLL2(theta, count, i, j, eLL);
                gradient(i,j) = negGrad;
                % gradient(i,j) = (c / theta(i,j)) - ((k - c) / (1 - theta(i,j)));
            end
        end
        R2 = R2 + gradient;
    end
end

%% Graph Alignment
kronExp = 4;
P = thetaG;
for i=1:kronExp
    P = kron(thetaG,P);
end
[p, pp] = firstPermutation(A, thetaG, 75000);

% Check graph alignment of A into P using the perms
n = size(A,1);
Ap = zeros(n,n);
for i=1:n
    for j=1:n
        Ap(i,j) = A(p(i), p(j));
    end
end

figure; 
subplot(1,3,1); imagesc(A); title('Kronecker Graph');
subplot(1,3,2); imagesc(P); title('Kronecker Expansion');
subplot(1,3,3); imagesc(Ap); title('Aligned Graph');

h = figure;
% h.Visible = 'off';
M(size(pp,1)) = struct('cdata',[],'colormap',[]);
for t = 1:size(pp,1)
    Ap = zeros(n,n);
    for i=1:n
        for j=1:n
            Ap(i,j) = A(pp(t,i), pp(t,j));
        end
    end
    imagesc(Ap);
    drawnow
    M(t) = getframe;
end


%%
clear all; clc; close all;
% filePath = "C:\Joshua/Software/snap/examples/as20graph.txt"
% filePath = "C:\Users\picka\Documents/software/snap/examples/as20graph.txt";
filePath = "C:\Joshua/Software/snap/examples\as20graph.txt";
E = readAdjList(filePath, 4);
A = sparse(E(:,1), E(:,2), 1);
theta0 = [0.9 0.6; 0.6 0.1];
theta0 = [0.895411 0.596948;
          0.596948 0.099491]
[theta, likelihoods] = NaiveKronFit(A, true, true, 2, theta0, 10);

%% Kronecker Hypergraph

HG = HAT.uniformErdosRenyi(27, 10, 3);
A = HG.adjTensor
A = (A > 0)
theta = rand(3,3,3);

HyperKronFit(A, v, true, 3, theta, 10)

inputs = repmat({[0,1]}, 1, k);

outputs = cell(1, k)
[outputs{:}] = ind2sub(size(A), linIdx)

(inputs{:})

[outputs{:}] = ndgrid(inputs{:})

%% 

eps = 1e-3;
theta = ones(2,2,2);
theta(2,2,2) = 0;
theta = theta - eps; theta(theta < 0) = eps;
A = superkron(theta, superkron(theta, theta))
linIdxs = 1:size(A,1)^ndims(A);
[x, y, z] = ind2sub(size(A), linIdxs');
data = [x, y, z, reshape(A,[numel(A), 1])]
figure;
scatter3(data(:,1),data(:,2),data(:,3),40,data(:,4),'filled')    % draw the scatter plot


%%

pts = 100000;
data = zeros(pts, 4);
for i=1:pts
    x = rand(2,1);
    y = rand(2,1);
    data(i,:) = kron(x,y);
end

figure;
labels = ["x_1","x_2","x_3","x_4"];
[h,ax] = plotmatrix(data);                        % create a 4 x 4 matrix of plots
for i = 1:4                                       % label the plots
  xlabel(ax(4,i), labels{i});
  ylabel(ax(i,1), labels{i});
end

figure;
scatter3(data(:,1),data(:,2),data(:,3),40,data(:,4),'filled')    % draw the scatter plot
xlabel("x_1"); ylabel("x_2"); zlabel("x_3"); title('Color is x_4')

%% Synthetic Test Graph 1
theta0 = [0.9 0.6;
          0.6 0.3];
n0 = 2;
kronExp = 14;
numE = 10 * n0^kronExp;
theta00 = theta0 / sum(theta0, 'all');

% Generate kronecker graph
E = kronGen(theta00, kronExp, numE);
writematrix(E, 'syntheticTestGraph4.txt', 'Delimiter', '\t')
% Graph 1 has kronExp = 8
% Graph 2 has kronExp = 10
% Graph 3 has kronExp = 12
% Graph 4 has kronExp = 14

n = n0^kronExp;
A0 = sparse(n, n);
for i=1:size(E,1)
    A0(E(i,1), E(i,2)) = 1;
end
A0 = full(A0);

figure; imagesc(A0)


