%% Calculating Kronecker Factorization of Toy Hypergraphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 11, 2023

clear; clc; close all;

n = 49;
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

% saveas(gcf, "Kronecker Factors of Toy Hypergraphs.png");