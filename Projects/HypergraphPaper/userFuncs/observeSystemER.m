function [O] = observeSystemER(HG,thresh,itrs,rH,rL)
%OBSERVESYSTEM 
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

r = [rH, 0, rL];
O = cell(itrs, 2*length(r));

for i=1:itrs
    for HOMB=1:6
        % Select hypergraph
        if i == 1
            HGa = HG;
        else
            HGa = O{i-1,HOMB};
        end
        [~, m] = size(observedHG(HGa).IM);
        mp = floor((1-thresh) * m);
%        disp(mp);
        % Perform observation on the system
        if HOMB == 1
            O{i,HOMB} = observeV(HGa,r(1),mp);
        elseif HOMB == 2
            O{i,HOMB} = observeV(HGa,r(2),mp);
        elseif HOMB == 3
            O{i,HOMB} = observeV(HGa,r(3),mp);
        elseif HOMB == 4
            O{i,HOMB} = observeE(HGa,r(1),mp);
        elseif HOMB == 5
            O{i,HOMB} = observeE(HGa,r(2),mp);
        elseif HOMB == 6
            O{i,HOMB} = observeE(HGa,r(3),mp);
        end
%        disp(string(HOMB) + ": " + string(size(observedHG(O{i,HOMB}).IM)));
    end
%    disp("======")
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
