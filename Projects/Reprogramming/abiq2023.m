function [outputArg1,outputArg2] = abiq2023(D)
%ABIQ2023 Sensor selection from Learning perturbation-inducible cell states
% from observability analysis of transcriptome dynamics
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 13, 2023

[n,t] = size(D);

%% Dynamic Mode Decomposition
out = DMD(D, [], 0.9);
Atilda = out.DMD.UX' * out.Xp * out.DMD.VX * inv(out.DMD.Sig);

%{
[e,v] = eig(Atilda)
figure;
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p, 'FaceColor', 'r'); hold on;
scatter(real(diag(v)), imag(diag(v)))
axis equal
%}
% plot(diag(v))
% norm(diag(v) - real(diag(v)))
% A = exactDMD(D);
% [e,v] = eig(A);

%% Reduced Grammarian

Ai = Atilda;
G = out.DMD.UX' * D(:,1) * D(:,1)' * out.DMD.UX;
for i=1:t
    G = G + Ai * out.DMD.UX' * D(:,i) * D(:,i)'*out.DMD.UX*Ai';
    Ai = Ai * Atilda;
end

[V,~] = eig(G);         % Eigen decomposition of reduced Grammarian
GV = out.DMD.UX * V;    % Approximate eigenvectors of full Grammarian

%% Sensor Rankings
[~, s] = max(GV);       % Take the largest index in each eigenvector

end

%{

    Gproj = out.DMD.UX * G * out.DMD.UX';
    
    [V,D] = eigs(Gproj,1)
    [m, mi] = max(V)
    
    issymmetric(Gproj)

%}
