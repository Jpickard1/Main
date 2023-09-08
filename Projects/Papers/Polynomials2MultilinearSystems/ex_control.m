%% Examples of Controllable Nonhomogeneous Systems
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 1, 2023

%% Example 2 of Polynomial Control Systems (this system is controllable)
p = mvpoly('type', 'JC2');
m = multilinearSystem('poly',p);
C = ctrbkc(double(m.A), 1);

% Controllability Test
if rank(C) == p.nvars
    disp("The system is strongly controllable")
end

S = spanOfLieAlg(p, [1; 0])
rank(S) == rank(C)

%% Example 2.2 of Controllability and Observability of Polynomial Dynamical
% Systems (this system is controllable)

p = mvpoly('type', 'JB22');
m = multilinearSystem('poly',p);
C = ctrbkc(double(m.A), 1);

% Controllability Test
if rank(C) == p.nvars
    disp("The system is strongly controllable")
else
    disp("The system is NOT strongly controllable");
end

S = spanOfLieAlg(p, [1; 0])
rank(S) == rank(C)

%% Example: randomly generated polynomial

n = 5;
d = 3;
controlNodes = 1;
esv = numPolyTerms(n, d);
Am = rand(n,esv);
p = mvpoly('Am',Am,'maxD',d);

m = multilinearSystem('poly',p);
C = ctrbkc(double(m.A), controlNodes);

S = spanOfLieAlg(p, getBmatrix(controlNodes,n))




