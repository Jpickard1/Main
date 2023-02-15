% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

N = 50; E = 100; K = 3;
rL = -1; rH = 1;
thresh = 0.05;
itrs = 5;

HG1 = HAT.uniformErdosRenyi(N,E,K);
O1 = observeSystemER(HG1,thresh,itrs,rH,rL);
HG2 = HAT.uniformErdosRenyi(N,E,K);
O2 = observeSystemER(HG2,thresh,itrs,rH,rL);
HG3 = HAT.uniformErdosRenyi(N,E,K);
O3 = observeSystemER(HG3,thresh,itrs,rH,rL);
HG4 = HAT.uniformErdosRenyi(N,E,K);
O4 = observeSystemER(HG4,thresh,itrs,rH,rL);
HG5 = HAT.uniformErdosRenyi(N,E,K);
O5 = observeSystemER(HG5,thresh,itrs,rH,rL);
O = [O1; O2; O3; O4; O5];
T = {HG1, HG2, HG3, HG4, HG5};

t = tableMaker(O', 'U-ER', itrs, 5, T, true,["OVH","OVN","OVL","OEH","OEN","OEL"]);

