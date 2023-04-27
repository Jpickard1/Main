% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 21, 2023

A = rand(5,5,5);
B = rand(5,5,5);
C = superkron(A, B);

A = tensor(A);
B = tensor(B);
C = tensor(C);

R = 2;
MA = cp_als(A,R);
MB = cp_als(B,R);
MC = cp_als(C,R^2);

MA.lambda
MB.lambda
(kron(MA.lambda, MB.lambda))'
MC.lambda'
