function glExpKronTenEig(n1, n2)
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 12, 2023

addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/'));
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/Dependencies/TenEig2/'));

itrs = 1000;
tB = zeros(itrs,1);
tC = zeros(itrs,1);
tA = zeros(itrs,1);
BCAerrors = zeros(itrs,1);
parfor i=1:itrs
    disp(i)
    B = rand(n1,n1,n1);
    C = rand(n2,n2,n2);
    A = superkron(B,C);
    
    tic;
    eb = heig(B);
    tB(i) = toc;
    tic;
    ec = heig(C);
    tC(i) = toc;
    tic;
    ea = heig(A);
    tA(i) = toc;

    BCAerrors(i) = max(ea) - max(eb)*max(ec);

end

save("eigencalculationRunTimes" + string(n1*n2) + ".mat")

end

