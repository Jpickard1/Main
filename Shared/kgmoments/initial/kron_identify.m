%% Identifiability of Kronecker graphs

%% Setup
addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl/');


%% Parameters


Ts = {
    [0.992 0.57 0.05],
};

ks = [15];
nrep = 4;

%%

for Ti=1:length(Ts)
    T = Ts{Ti};
    a = T(1);
    b = T(2);
    c = T(3);
    for ki=1:length(ks)
        k = ks(ki);
        name = sprintf('kron-%3.1f-%3.1f-%3.1f-%i',a,b,c,k);
        for ri=1:nrep
            G = rmat(k,[a b; b c]); 
            % symmetrize the graph
            G = triu(G,1);
            G = G|G';
            props = graph_properties(G,name);
            [params,obj,G]= kron_moment_fit(props);
        end
    end
end
           

    