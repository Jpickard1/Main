clear; clc;

n = 250;
P = rand(n,n);
A = zeros(size(P));
for i=1:n
    for j=1:n
        A(i,j) = binornd(1,P(i,j));
    end
end

figure;
scatter(1:n, sort(real(eig(P)), 'descend'), '.'); hold on;
scatter(1:n, sort(real(eig(A)), 'descend'), '.');

max(eig(P))
max(eig(A))
