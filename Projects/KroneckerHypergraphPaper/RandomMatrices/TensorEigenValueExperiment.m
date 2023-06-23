function TensorEigenValueExperiment(sym, worker)
%TENSOREIGENVALUEEXPERIMENT
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 23, 2023

n = 4;
T = zeros(n,n,n);

for i=1:n
    for j=1:n
        for k=1:n
            element = rand();
            T(i,j,k) = element;
            T(i,k,j) = element;
            T(j,i,k) = element;
            T(j,k,i) = element;
            T(k,i,j) = element;
            T(k,j,i) = element;
        end
    end
end

heigFull = heig(T)

n = 2^10; e = 2*n; k = 2;
HG = HAT.uniformErdosRenyi(n, e, k)
E = HAT.uniformEdgeSet(HG);
A = HG.adjTensor;
theta0 = rand(2,2)

[theta] = HyperKronFit('E', E, 'v', true, 'gradSamples', 100000, 'firstPermItrs',1000, 'learningRate',1e-5, 'maxItrs',30, 'theta0', theta0)


%% Time large heig calculation
n = 100;
A = zeros(n,n,n);
for i=1:n
    for j=i:n
        for k=j:n
            element = rand();
            A(i,j,k) = element;
            A(i,k,j) = element;
            A(j,i,k) = element;
            A(j,k,i) = element;
            A(k,i,j) = element;
            A(k,j,i) = element;            
        end
    end
end
tic;
A = heig(A,'symmetric')
t = toc;
disp(toc);

end

