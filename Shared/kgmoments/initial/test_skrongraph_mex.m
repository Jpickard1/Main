T = [0.5 0.3; 0.3 0.1];
P = kron(T,T);
P = kron(T,P);

n = 8;
ntrials = 1000000;
As = zeros(n,n);
for t = 1:ntrials
    edges = skrongraph_mex(T,3,0);
    A = sparse(edges(:,1),edges(:,2),1,n,n);
    As = As + A;
end

As./ntrials
P