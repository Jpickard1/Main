%% Identifiability of Kronecker graphs

%% Setup
addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl/');
addpath('../../matlab');


%% Parameters


Ts = {
    [0.99,0.48,0.25], % 
    [1.00,0.67,0.08],
    2.32*[0.45,0.15,0.25],
    2.32*[0.57,0.19,0.05],
    [0.992,0.57,0.05],
    [0.999,0.271, 0.587],
};


ks = [14];
nrep = 10;

%%

for Ti=1:length(Ts)
    T = Ts{Ti}
    a = T(1);
    b = T(2);
    c = T(3);
    for ki=1:length(ks)
        k = ks(ki);
        name = sprintf('kron-%3.1f-%3.1f-%3.1f-%i',a,b,c,k);
        for ri=1:nrep
            G = rmat([a b; b c],k); 
            % symmetrize the graph
            G = triu(G,1);
            G = G|G';
            props = graph_properties(G,name);
            [params,obj,G]= kron_moment_fit(props);
        end
    end
end
    
