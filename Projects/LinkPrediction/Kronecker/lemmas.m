%% Lemmas
%
%   In this file I numerically verity several lemmas I am interested in for
%   my Kronecker link prediction paper.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 2, 2023

%% Lemma 1
% The degree sequence of a matrix can be written as the degree sequene of
% its Kronecker factors

% Set dimensions and create factor matrices a, b, and A
n1 = 2;
n2 = 3;
m1 = 2;
m2 = 3;
a  = rand(n1, m1);
b  = rand(n2, m2);
A = kron(a, b);

% Degree sequence of A
dA = sum(A, 2);

% Get degree sequence of A from a and b
Dab = [];
for i=1:size(Aab, 1)
    Dab(i) = sum(kron(a(ceil(i/n2), :), b(mod(i-1, n2)+1,:)));
end

dA'
Dab

