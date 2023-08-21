%% SIS System
%
%   This file constructs the SIS system as a multilinear system.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: August 17, 2023

%% Classic SIS System (mvpoly)
beta = 0.25;
gamma = 0.5;
p = mvpoly('type','SIS','beta',beta,'gamma',gamma,'stoch',0);
[t,X] = ode45(@(t, x) p.eval(x),[0 100],[0.25 0.75]);
figure; hold on;
plot(t,X(:,1));
plot(t,X(:,2));
title(p.title(),'Interpreter','latex');

%% Distribution of Eigenvalues for a stochastic process

beta = 0.25;
gamma = 0.5;
itrs = 50;
Emax = zeros(1,itrs);
R = 0.1*rand(itrs,1);

for i=1:itrs
    disp(string(i) + ": " + string(R(i)));
    p = mvpoly('type','SIS','beta',beta,'gamma',gamma,'stoch',R(i));
    m = multilinearSystem('poly',p);
    A = double(m.A);
    Emax(i) = max(zeig(A));
end

figure; histogram(100*Emax)
figure; plot(sort(Emax))

