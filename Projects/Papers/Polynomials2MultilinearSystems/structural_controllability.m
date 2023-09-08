clear all; close all; clc;

A = zeros(2,2,2);
A(1,1,1) = 1;
A(2,2,1) = 1;
A(1,1,2) = 1;
A(2,2,2) = 1;

A

HAT.ctrbk(A, 1)

%%
clear all; close all; clc;

A = zeros(2,2,2);
A(1,1,1) = 1;
A(2,1,1) = 1;
A(1,1,2) = 1;
A(2,2,2) = 1;

A

reshape(A, [2 4])

rank(HAT.ctrbk(A, 1))

