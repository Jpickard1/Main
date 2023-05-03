%% Kronecker Hypergraph Adjacency Tensor
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 3, 2023

H = hyperring(5, 3);
A1 = H.adjTensor;

idxs = find(A1 ~= 0);
[x, y, z] = ind2sub(size(A1), find(A1 ~= 0));
% figure; scatter3(x,y,z);

A = A1;
itrs = 4;
figure;
for itr=1:itrs
    if itr == 1; A = A1; else; A = superkron(A, A1); end;
    [x, y, z] = ind2sub(size(A), find(A ~= 0));
    subplot(1, itrs, itrs - itr + 1); scatter3(x,y,z, '.');
end
