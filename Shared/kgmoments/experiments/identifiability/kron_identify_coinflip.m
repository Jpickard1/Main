%% Identifiability of Kronecker graphs
% In this experiment

%% Setup
addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl/');
addpath('../../matlab');


%% Parameters



Ts = {
    [0.99,0.48,0.25] % 
    [1.00,0.67,0.08]
    %[0.992,0.57,0.05],
    [0.999,0.271, 0.587]
    [0.87,0.6,0.7]
    2.32*[0.57,0.19,0.05]
};

ks = [14];
nrep = 50;


%%
rand('state',1);
randn('state',1);
results = [];
for Ti=1:length(Ts)
    T = Ts{Ti};
    a = T(1);
    b = T(2);
    c = T(3);
    for ki=1:length(ks)
        k = ks(ki);
        name = sprintf('kron-%5.3f-%5.3f-%5.3f-%i',a,b,c,k);
        for ri=1:nrep
            tic;
            G = skrongraph([a b; b c],k,'nthreads',8); 
            % symmetrize the graph
            G = triu(G,1);
            G = G|G';
            props = graph_properties(G,name);
            [params,stats,G]= kron_moment_fit(props,'sampletype','skron','obj','var');
            toc;
            result.props = props;
            result.T = T;
            result.Ti = Ti;
            result.a = a;
            result.b = b;
            result.c = c;
            result.k = k;
            result.ki = ki;
            result.rep = ri;
            result.params = params;
            
            result.stats = stats;
            result.obj = stats.obj;
            
            if ~isempty(results)
                results(end+1) = result;
            else
                results = result;
            end
            save 'coinflip_identify.mat' results Ts nrep ks;
        end
    end
end
           

    