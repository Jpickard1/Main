%% SIR Models: How can the Kronecker framework be applied here?
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 26, 2023

B = [1 1;
     1 0];
C = [0 1 0;
     1 0 1;
     0 1 0];

A = kron(B,C);

figure;
subplot(1,3,1); plot(graph(A)); title('A');
subplot(1,3,2); plot(graph(B)); title('B');
subplot(1,3,3); plot(graph(C)); title('C');

