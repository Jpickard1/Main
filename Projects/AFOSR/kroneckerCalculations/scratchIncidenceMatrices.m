%% Incidence Matrices
%
%   How can the incidence matrix of a Kronecker hypergraph be expressed in
%   terms of its factors?
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 24, 2023

clc; clear all; close all;

IM = [1 1 1 0;
      0 1 1 1];
IM = IM';

H1 = Hypergraph('IM', IM);
A1 = H1.adjTensor;
V1 = ["x_1","x_2","x_3","x_4"];
V1y = ["y_1","y_2","y_3","y_4"];

Ak = superkron(A1, A1);
IMk = HAT.A32IM(Ak);
Hk = Hypergraph('IM', IMk);
Vk = kronString(V1, V1y);

% figure;
% subplot(1,2,1); H1.plot(); yticks(1:4); yticklabels(V1);
% subplot(1,2,2); Hk.plot(); yticks(1:16); yticklabels(Vk);

figure;
subplot(1,2,1); HAT.plotIncidenceMatrix(H1, 'sort', false); yticks(1:4); yticklabels(V1);
subplot(1,2,2); HAT.plotIncidenceMatrix(Hk, 'sort', false); yticks(1:16); yticklabels(Vk);

for e=1:12 %size(IMk,2)
    vx = find(IMk(:,e) ~= 0);
    str = "";
    for i=1:length(vx)
        str = str + "(" + Vk(vx(i)) + ") "; 
    end
    disp(str);
end

% >>> (x_1 - y_1) (x_2 - y_2) (x_3 - y_3) 
% >>> (x_1 - y_1) (x_2 - y_3) (x_3 - y_2) 
% >>> (x_1 - y_2) (x_2 - y_1) (x_3 - y_3) 
% >>> (x_1 - y_2) (x_2 - y_3) (x_3 - y_1) 
    % >>> (x_1 - y_2) (x_2 - y_3) (x_3 - y_4) 
    % >>> (x_1 - y_2) (x_2 - y_4) (x_3 - y_3) 
% >>> (x_1 - y_3) (x_2 - y_1) (x_3 - y_2) 
% >>> (x_1 - y_3) (x_2 - y_2) (x_3 - y_1) 
    % >>> (x_1 - y_3) (x_2 - y_2) (x_3 - y_4) 
    % >>> (x_1 - y_3) (x_2 - y_4) (x_3 - y_2) 
    % >>> (x_1 - y_4) (x_2 - y_2) (x_3 - y_3) 
    % >>> (x_1 - y_4) (x_2 - y_3) (x_3 - y_2) 
% 
% OUTPUT: is all possible combinations of (x_1, x_2, x_3) with (y_1, y_2,
% y_3) and (y_2, y_3, y_4) of which there are 12. 
%
% How do we understant this? Well, there are 6 new hyperedges composed of
% (x_1, x_2, x_3) and one hyperedge from y (I split these up in the output 
% with indentations). To create each new hyperedges, we must select each
% possible pairing of the x and y vertices. With the x vertices in a fixed
% order, the number of possible pairings is equivalent to the number of
% permutations of the y vertices, and for 3-uniform hyperedges, there are 6
% possible orderings.
%
% i.e. the numeber of hyperedges in HGK = (number of hyperedges in factor
% 1) x (number of hyperedges in factor 2) x (size of permutaiton set of a 
% single hyperedge), where factor 1 and factor 2 are both k uniform
% hypergraphs.
%
%% How can I write this in some matrix form?
%
%   I could generate all possible permutations of a hyperedge by stacking
%   all possible permutation matrices.

n=3;
E = eye(n);
pms = perms(1:n);
P1 = E(pms(1,:),:);
P2 = E(pms(2,:),:);
P3 = E(pms(3,:),:);
P4 = E(pms(4,:),:);
P5 = E(pms(5,:),:);
P6 = E(pms(6,:),:);
X = sym("x_%d", [n 1]);
XP = reshape([P1 P2 P3 P4 P5 P6]' * X, [3 6])

%% What happens when we take the Kronecker product of the IMs?

kron(IM, IM)
