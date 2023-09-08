%% Homogeneous Polynomial Vector 
%
% Jurdjevic and Kupka: a polynomial field is called homogeneous of degree P
% if P(lambda x) = lambda^p P(x)

x = 2;
y = 2;
lambda = 3;

n = 3;
d = 2;
esv = numPolyTerms(n, d);
Am = sparse(n, esv);
%     xx xy xz yx yy yz zx zy zz  x  y  z  c
%      1  2  3  4  5  6  7  8  9 10 11 12 13
Am(1,1) = rand();
Am(1,2) = rand();
Am(1,3) = rand();
Am(2,4) = rand();
Am(2,5) = rand();
Am(2,6) = rand();
Am(3,7) = rand();
Am(3,8) = rand();
Am(3,9) = rand();

p = mvpoly("Am",Am,'nvars',n,'maxD',d)
x = rand(3,1);
lambda = rand();
p.eval(lambda * x)
lambda^d * p.eval(x)

%%
sigma = 10;
rho = 28;
beta = 8/3;
p = mvpoly('type','lorenz','sigma',sigma,'rho',rho,'beta',beta);
m = multilinearSystem('poly',p);

x = rand(3,1);
lambda = rand();
m.evalLambda(x, lambda)
lambda^2 * m.eval(x)


