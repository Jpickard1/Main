n = 20;
B = rand(20,20,20);
A = superkron(B,B);
tic;
aa = tt_tensor(B);
disp(toc)
tic;
aa = tt_tensor(A);
disp(toc)

n = 20;
A = rand(n,n,n);
A = tensor(A);
tol = 1e-4;
r = 4;
tic;
cp_als(A, r);
tt1 = toc
A = rand(n,n,n);
A = superkron(A,A)
A = tensor(A);
tic;
cp_als(A, r);
tt2 = toc
