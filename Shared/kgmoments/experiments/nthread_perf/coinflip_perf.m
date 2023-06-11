%% Multithreaded performance of the coinflipping procedure
% A precise way to generate stochastic kronecker graphs is to use explicit
% coinflipping.  This avoids the approximation in the number of edges.  
% In this experiment, we investigate how many threads to use.  

%% Results summary
% This summary is based on the experiment on David's Core i7-960 on April
% 10th, 2011 on Ubuntu 10.04, Matlab 2010b.
% 
% It shows that we continue to get results faster using ALL cores.  Wild!
% So hyperthreading actually helps here.

%% The experiment
% Kronecker parameters
% these don't really matter.
a = 0.992;
b = 0.57;
c = 0.05;

% the experimental conditions
% These do matter.
nthreads = 1:8;
nrep = 5;
k = 12;

for t=nthreads
    fprintf('%i thread\n',t);
    tic, for i=1:5, tes = skrongraph_mex([a b; b c],k,t); end, toc
end