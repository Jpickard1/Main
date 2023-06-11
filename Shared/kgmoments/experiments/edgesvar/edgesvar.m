%% Check the variance in the number of edges in a coinflip model.

%% Setup
addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl/');
addpath('../../matlab');



%%
Ts = {
    [0.99,0.48,0.25], % 
    [1.00,0.67,0.08],
    2.32*[0.45,0.15,0.25],
    2.32*[0.57,0.19,0.05],
    [0.992,0.57,0.05],
    [0.999,0.271, 0.587],
};


k = 12;
nrep = 50;
%%
nedges = zeros(length(Ts),nrep);
for Ti=1:length(Ts)
    T = Ts{Ti}
    a = T(1);
    b = T(2);
    c = T(3);
  
    for ri=1:nrep
        G = skrongraph([a b; b c],k,'nthreads',8); 
        % symmetrize the graph
        G = triu(G,1);
        G = G|G';
        nedges(Ti,ri) = nnz(G)/2;
    end
end

%%
hist(nedges(1,:),10)
%%
hist(nedges(2,:),10)
%%
hist(nedges(3,:),10)
%%
hist(nedges(4,:),10)
%%
hist(nedges(5,:),10)

%% 
% These results show that the result is roughly a normal distribution.  The
% variance is about the square-root of the expectation.