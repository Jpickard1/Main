clear all; close all; clc

tol = 1e-6;
n1 = 10;
n2 = 10;

B = rand(n1,n1,n1);
C = rand(n2,n2,n2);
A = superkron(B, C);

B = tensor(B);
C = tensor(C);
A = tensor(A);

Bh = hosvd(B, tol);
Ch = hosvd(C, tol);
Ah = hosvd(A, tol);

BhC = double(Bh.core);
ChC = double(Ch.core);
AhC = double(Ah.core);

BhCKChC = superkron(BhC, ChC);

sum(abs(AhC) - abs(BhCKChC), 'all') % This line
sum(AhC - BhCKChC, 'all') % This line
sum((AhC), 'all')

%%

n1 = 5; n2 = 5;
A = rand(n1,n1,n1);
B = rand(n2,n2,n2);
C = rand(n1,n1,n1);
D = rand(n2,n2,n2);

ACkBD = superkron(tensorprod(A,C), tensorprod(B,D));
AkBCkD = tensorprod(superkron(A,B), superkron(C,D));

max(abs(ACkBD - AkBCkD), [], 'all')


