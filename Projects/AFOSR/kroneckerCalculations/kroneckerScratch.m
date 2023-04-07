% Auth: Joshua Pickard
%       jpic@umich.edu

%% What is bigger: log(x) + log(y) or log(x)log(y) for x, y < 0

itrs = 1000;
c = 0;
for i=1:itrs
    x = rand(); y = rand();
    if log(x*y) > - log(x) * log(y)
        c = c + 1;
    end
end
disp(c)

%% What does the unfolding of Kronecker tensors look like?
%   It is not clear to me yet
clear
x = sym("x",[2 2 2]);
y = sym("y",[2 2 2]);

x = rand([2,2,2]);
y = rand([2,2,2]);
x1 = reshape(x, [2 4]);
y1 = reshape(y, [2 4]);
x1y1 = kron(x1,y1)
xy = superkron(x,y);
xy1 = reshape(xy, [4 16])

xy1 == x1y1

%% Is is true that kron(x^[m-1], y^[m-1]) = kron(x,y)^[m-1]
% >> Yes
for i=1:100
n1 = 5;
n2 = 5;
m = 3;
x=rand(n1,1); y=rand(n2,1);
x2 = x .^ 2; y2 = y .^ 2;
x2y2 = kron(x2,y2);
xy = kron(x,y); xy2 = xy .^ 2;

if sum((xy2 - x2y2) > 1e-5) > 0
    disp('Learn')
end
% disp(i)
end

%% Date: April 5, 2023
close all
HG1 = getToyHG(3, 3, 'hyperchain'); A1 = HG1.adjTensor;
HG2 = getToyHG(3, 3, 'hyperchain'); A2 = HG2.adjTensor;

A = superkron(A1, A2);

idxs = find(A > 0);
[i1 i2 i3] = ind2sub(size(A), idxs)
E = [i1 i2 i3]; E = sortrows(E);
E = unique(E,'rows');
IM = HAT.hyperedge2IM(E);
HG = Hypergraph('IM', IM);
figure; HG.plot()
figure; HG1.plot()
figure; HG2.plot()
