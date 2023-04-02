%% Complete Hypergraphs
%
%   Here I compute the MON on 3-uniform complete hypergraph with 5 and 6
%   vertices respectively
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 30, 2023

%% 3-uniform complete hypergraph
n = 5;
k = 3;
ES = nchoosek(1:n, 3);
IM = HAT.hyperedge2IM(ES);
HG = Hypergraph('IM', IM);
O = HGObsvSym(HG);
[D, ~] = greedyMON(O, n);             % Greedy Node Selection
disp(D)
%>> disp(D)
%     1

%%
n = 6;
k = 3;
ES = nchoosek(1:n, 3);
IM = HAT.hyperedge2IM(ES);
HG = Hypergraph('IM', IM);
O = HGObsvSym(HG);
[D, ~] = greedyMON(O, n);             % Greedy Node Selection
disp(D)
%>> disp(D)
%   1

%% figure out recursive size of sets
k = 3;
P = 10;
p = P;
S = (p-1)*k-(2*p-3);
for p=9:-1:1
    S = S * (p-1)*k-(2*p-3);
    disp(S)
end

