%% 03312023 - Revieiwng results used in IEEE CDC 2023 Paper submission
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 31, 2023

%% Quick claculations on hyperstar(5,3)
HG = getToyHG(5,3,'hyperstar');
O = HGObsvSym(HG)
disp(O{1})
a = O{1};
rank(O{3})
rank([O{3}; O{1}]) % >> 4
rank([O{3}; O{4}]) % >> 5
rank(O{3} - O{1}) % >> 4
rank(O{3} - O{4}) % >> 1

%% Check recursive size for Jp3
k = 3;
r = (6-1)*3-(12-3);
for p=5:-1:2
    b = (p-1)*k-(2*p-3);
    r = r * b; disp(r)
end

%% testing Jp3
clear; clc;

n = 5; k = 3;
HG = getToyHG(n, k, 'hyperstar');


p = n;
rpts = p*k-(2*p-1);
loops = (p-1) * k- (2 * p - 3);
S = cell(loops, 1);

for i=1:loops
    S{i} = repmat(xx, 1, rpts);
end
rpts = p*k-(2*p-1);
xx = sym('x', [n, 1]);
S = repmat(xx, 1, rpts);
Jp3(HG, n, S)

%% Recalculating results with Jp3
n = 5; k=3;
HG = getToyHG(n, k,'hyperstar');
O1 = HGObsvSym(HG)
O0 = HGObsvSym0(HG)
[D, ~] = greedyMON(O, n);             % Greedy Node Selection


O{1}(2,:)

OO = [];
for i=1:5
    disp(rank(O{i}));
    OO = [OO; O{i}];
end

rank([O{4}; O{3}])

O{1}(1,:)

rank(OO)

