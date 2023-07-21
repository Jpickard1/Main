function glExpEig()
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 15, 2023

addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/'));
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/Dependencies/TenEig2/'));

itrs = 5;
maxN = 20;
tB = zeros(maxN,itrs);
tC = zeros(maxN,itrs);
tA = zeros(maxN,itrs);
cB = cell(maxN,itrs);
cC = cell(maxN,itrs);
cA = cell(maxN,itrs);

BCAerrors = zeros(itrs,1);
for j=1:20
    disp(j);
    parfor i=1:itrs
        disp(i)
        B = rand(n1,n1,n1);
        C = rand(n2,n2,n2);
        A = superkron(B,C);
        
        tic;
        eb = heig(B);
        tB(j,i) = toc;
        tic;
        ec = heig(C);
        tC(j,i) = toc;
        tic;
        ea = heig(A);
        tA(j,i) = toc;

        cB{j,i} = eB;
        cC{j,i} = eC;
        cA{j,i} = eA;
    end
    save("heigValues.mat")
end
end

