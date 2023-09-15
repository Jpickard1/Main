function [A] = SDMD(DataMatrix, S)
%SDMD Structured Dynamic Mode Decomposition
%
%   min_A |Xp - X * A| s.t. (A ~= 0) == S
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 14, 2023

% Construct data matrices note, time snapshots are columns
    % X  = DataMatrix(:,1:end-1);
    % Xp = DataMatrix(:,2:end);

% Compute Xp' \ X' 
%   min |Xp' - X'*A| s.t. (A ~= 0) == S
%     A = structuredLS(Xp',X',S');
A = zeros(size(S));
C = zeros(size(S));
for i=1:size(C,2)
    N = find(S(:,i) ~= 0);
    N = N(N ~= i);
    Ni = [i; N];
    out = DMD(DataMatrix(Ni,:),[],0.9);
    A_bar = out.DMD.A_bar;
    A(N,i) = A_bar(1,2:end);
    % Ci = B(:,N) \ A(:,i);
    % C(N,i) = Ci;
end

end

