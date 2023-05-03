%% Calculating Kronecker Factorization of Toy Hypergraphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 11, 2023

clear; clc; close all;

n = 25;
k = 3;

figure;
HG = hyperchain(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,1); HG.plot(); title('Hyperchain');
        subplot(3,3,2); HG1.plot(); title('Factor 1');
        subplot(3,3,3); HG2.plot(); title('Factor 2');
HG = hyperring(n,k);
A = HG.adjTensor;   % A = tensor(A);
[B, C] = NTKP(A) ;  %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,4); HG.plot(); title('Hyperring');
        subplot(3,3,5); HG1.plot(); title('Factor 1');
        subplot(3,3,6); HG2.plot(); title('Factor 2');
HG = hyperstar(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,7); HG.plot(); title('Hyperstar');
        subplot(3,3,8); HG1.plot(); title('Factor 1');
        subplot(3,3,9); HG2.plot(); title('Factor 2');

saveas(gcf, "Kronecker Factors of Toy Hypergraphs " + string(n) + ".png");

%% Adjacency plots
figure;
HG = hyperchain(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,1);  [x, y, z] = ind2sub(size(A), find(A ~= 0)); scatter3(x,y,z, '.'); title('Hyperchain'); hold on;
        subplot(3,3,2); [x, y, z] = ind2sub(size(B), find(B ~= 0)); scatter3(x,y,z, '.');  title('Factor 1'); hold on; 
        subplot(3,3,3); [x, y, z] = ind2sub(size(C), find(C ~= 0)); scatter3(x,y,z, '.');  title('Factor 2'); hold on; 
HG = hyperring(n,k);
A = HG.adjTensor;   % A = tensor(A);
[B, C] = NTKP(A) ;  %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,4); [x, y, z] = ind2sub(size(A), find(A ~= 0)); scatter3(x,y,z, '.'); title('Hyperchain'); hold on;
        subplot(3,3,5); [x, y, z] = ind2sub(size(B), find(B ~= 0)); scatter3(x,y,z, '.');  title('Factor 1'); hold on; 
        subplot(3,3,6); [x, y, z] = ind2sub(size(C), find(C ~= 0)); scatter3(x,y,z, '.');  title('Factor 2'); hold on; 
HG = hyperstar(n,k);
A = HG.adjTensor; % A = tensor(A);
[B, C] = NTKP(A) ; %, [10 10 10], [10 10 10]);
HG1 = Hypergraph('IM', HAT.A32IM(B));
HG2 = Hypergraph('IM', HAT.A32IM(C));
        subplot(3,3,7); [x, y, z] = ind2sub(size(A), find(A ~= 0)); scatter3(x,y,z, '.'); title('Hyperchain'); hold on;
        subplot(3,3,8); [x, y, z] = ind2sub(size(B), find(B ~= 0)); scatter3(x,y,z, '.');  title('Factor 1'); hold on; 
        subplot(3,3,9); [x, y, z] = ind2sub(size(C), find(C ~= 0)); scatter3(x,y,z, '.');  title('Factor 2'); hold on; 
