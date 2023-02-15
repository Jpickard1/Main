function [O] = observeSystemLP(HG,thresh,r)
%OBSERVESYSTEM 
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 13, 2023

% Select hypergraph
HGa = HG;

% Get number vx and edges
[v, m] = size(observedHG(HGa).IM);
mp = floor((1-thresh) * m); % * v);

O = cell(2*length(r));

for HOMB=1:(length(O)/2)
    O{HOMB} = observeV(HGa,r(HOMB),mp);
    O{HOMB + length(O)/2} = observeE(HGa,r(HOMB),mp);
end

end

%% Testing
%{
N = 10; E = 100; K = 5;
rL = -1; rH = 1;
thresh = 0.03;
itrs = 20;
HG1 = HAT.uniformErdosRenyi(N,E,K)
O1 = observeSystemER(HG1,thresh,itrs,rH,rL);
HG2 = HAT.uniformErdosRenyi(N,E,K)
O2 = observeSystemER(HG2,thresh,itrs,rH,rL);
HG3 = HAT.uniformErdosRenyi(N,E,K)
O3 = observeSystemER(HG3,thresh,itrs,rH,rL);
HG4 = HAT.uniformErdosRenyi(N,E,K)
O4 = observeSystemER(HG4,thresh,itrs,rH,rL);
HG5 = HAT.uniformErdosRenyi(N,E,K)
O5 = observeSystemER(HG5,thresh,itrs,rH,rL);
O = [O1; O2; O3; O4; O5];
T = {HG1, HG2, HG3, HG4, HG5};

tableMaker(O', 'U-ER', itrs, 5, T, false,0);
%}

%{
% HG = HAT.uniformErdosRenyi(N,E,K)
% A = HG.adjTensor;
% hosvd(A, 1e-2,'Rank',[5,5,5])
IM = full(HG.IM);
%}

%{
    if HOMB == 1
%        mp = floor((1-thresh) * m * v);
        O{HOMB} = observeV(HGa,r(1),mp);
    elseif HOMB == 2
%        mp = floor((1-thresh) * m * v);
        O{HOMB} = observeV(HGa,r(2),mp);
    elseif HOMB == 3
%        mp = floor((1-thresh) * m * v);
        O{HOMB} = observeV(HGa,r(3),mp);
    elseif HOMB == 4
        O{HOMB} = observeE(HGa,r(1),mp);
    elseif HOMB == 5
        O{HOMB} = observeE(HGa,r(2),mp);
    elseif HOMB == 6
        O{HOMB} = observeE(HGa,r(3),mp);
    end
end
%}
