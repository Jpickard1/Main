
%% One thread
T = [0.5 0.3; 0.3 0.1];
P = kron(T,T);
P = kron(T,P);

n = 8;
ntrials = 1000000;
As = zeros(n,n);
for t = 1:ntrials
    %edges = skrongraph(T,3,1,);
    A = sparse(edges(:,1),edges(:,2),1,n,n);
    As = As + A;
end

As./ntrials
P


%% Two threads
T = [0.5 0.3; 0.3 0.1];
P = kron(T,T);
P = kron(T,P);

n = 8;
ntrials = 1000000;
As = zeros(n,n);
for t = 1:ntrials
    edges = skrongraph_mex(T,3,2);
    A = sparse(edges(:,1),edges(:,2),1,n,n);
    As = As + full(A);
end

As./ntrials
P

%% Four threads
T = [0.5 0.3; 0.3 0.1];
P = kron(T,T);
P = kron(T,P);

n = 8;
ntrials = 1000000;
As = zeros(n,n);
for t = 1:ntrials
    edges = skrongraph_mex(T,3,4);
    A = sparse(edges(:,1),edges(:,2),1,n,n);
    As = As + full(A);
end

As./ntrials
P

%% Seven threads
T = [0.5 0.3; 0.3 0.1];
P = kron(T,T);
P = kron(T,P);

n = 8;
ntrials = 1000000;
As = zeros(n,n);
for t = 1:ntrials
    edges = skrongraph_mex(T,3,7);
    A = sparse(edges(:,1),edges(:,2),1,n,n);
    As = As + full(A);
end

As./ntrials
P