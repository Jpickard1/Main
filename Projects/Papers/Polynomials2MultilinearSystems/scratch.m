%% Polynomials to Multilinear Systems
%
%   POLYNOMIAL REPRESENTATION: Amat * (x kron ... kron x)
%
%   MULTILINEAR REPRESENTATION: sptensor class by Kolda and Bader in tensor
%   toolbox
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 5, 2023

%% mvpoly construction
clear; close all; clc;

Am = rand(3,13);
p = mvpoly("Am",Am,"maxD",2);
x = rand(3,1);
sv = [kron(x,x);x;1];

disp(Am * sv)
disp(p.eval(x))

%% multilinear system construction from a polynomial
clear; close all; clc;

n = 2;
Am = rand(2,3);
p = mvpoly("Am",Am,"maxD",1);
m = multilinearSystem("poly",p);

x = rand(n,1);

p.eval(x)
m.eval(x)

%% Test Evaluation of multlinear system vs a polynomial
clear; close all; clc;

n = randi([2 5], 1);
d = randi([2 5], 1);
npt = numPolyTerms(n,d);
Am = rand(n, npt, 0.1);
p = mvpoly("Am",Am,"maxD",d);
m = multilinearSystem("poly",p);

x = rand(n,1);

yp = p.eval(x);
ym = m.eval(x);

assert(sum(abs(yp - ym)) < 1e-10)

%% Kronecker indices

% Example data
x = [1; 2; 3];
d = 4;
i = 7;
n = 3;
% Call the function to get the indices j_1, j_2, ..., j_d
indices = kroneckerIndices(n, d, 80);

disp(indices)

x = sym('x_%d',[3,1]);
KroneckerPower(x,4)

%% Subtensors of sparse tensors

% Create the 3x3x3x3 sparse tensor A
sizeA = [3, 3, 3, 3];
A = sptensor(sizeA);

% Create the 3x3x3 sparse tensor B
sizeB = [3, 3, 3];
B = sptensor(sizeB);
B(1,1,1) = 1;
B(2,2,2) = 2;
B(3,3,3) = 3;

% Set B as a sub-tensor of A at the first dimension (A(:,:,:,1) = B)
A(:,:,:,1) = B;


