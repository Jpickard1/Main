%% This file utilizes data and methods from E. Yeung
%
% Data files were obtained from Cooper
%   - time_series_rna_fold_changes.csv
%   - time_series_rna_Z_scores.csv
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 3, 2023

%% Read in data
% Data matrices D1,D2, and D3 are replicate experiments with 15 time points
% each. The i,j-th entry the the value of RNA for gene i at time j computed
% according to the fold changes.

T = readtable("time_series_rna_fold_changes.csv");
D1 = zeros(15,size(T,2)-1);
D2 = zeros(15,size(T,2)-1);
D3 = zeros(15,size(T,2)-1);
for t=1:15
    D1(t,:) = T{3*(t-1) + 1, 2:size(T,2)};
    D2(t,:) = T{3*(t-1) + 2, 2:size(T,2)};
    D3(t,:) = T{3*(t-1) + 3, 2:size(T,2)};
end

D1 = D1'; D2 = D2'; D3 = D3';

%% Exact DMD
A1 = exactDMD(D1);
A2 = exactDMD(D2);
A3 = exactDMD(D3);

%% Kalman Filter

C1 = zeros(1,size(A1,1)); C1(1) = 1;

Ts = -1;
sys = ss(A1,[],C1,[],Ts,'OutputName','y');  % Plant dynamics and additive input noise w
[kalmf,L,~,Mx,Z] = kalman(sys,[],[]);


%% Sensor Selection


%% Apply DMD

O1 = DMD(D1, [], 0.9)
O1.DMD
