%% Evaluate objective functions for fitting Kronecker parameters
% There are a few ways to formulate the moment matching objective.  We
% investigate these variations here.

%% Setup

addpath('../../matlab');
addpath('~/dev/matlab-bgl');

%% Load precomputed graph feature data
% This experiment uses precomputed graph feature data exclusively.
% To add graphs to the experiment, add graphs to the
% 'data/graph_features.m' script.  

load('../../data/graphs_props.mat');

%% Add one more graph using Kronecker parameters
T=[0.99,0.48,0.25];
a = T(1);
b = T(2);
c = T(3);
k=14;
rand('state',0);
randn('state',0);
G = skrongraph([a b; b c],k,'nthreads',8); 
G = triu(G,1);
G = G|G';
name = sprintf('kron-%5.3f-%5.3f-%5.3f-%i',a,b,c,k);
props = graph_properties(G,name);
results(name) = props;
%%
% The output of loading this script is a containers.Map type called
% results. We will just look at a subset of the graphs.  
graphs = {name,'ca-GrQc','as20000102',};

%% Iterative over graphs
% For each graph, we want to compute a fit using both approaches: direct
% minimization and grid search.  

objresults = containers.Map;

for gi=1:length(graphs)
    graph = graphs{gi};
    gdata = results(graph);
    
    fits = [];
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'obj','diff1');
    
    fits.diff1.params = params;
    fits.diff1.stats = stats;

    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'obj','var');
    
    fits.var.params = params;
    fits.var.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'obj','var2');
    
    fits.var2.params = params;
    fits.var2.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'obj','diff1var');
    
    fits.diff1e.params = params;
    fits.diff1e.stats = stats;

    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'obj','varf');
    
    fits.varf.params = params;
    fits.varf.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'obj','var2f');
    
    fits.var2f.params = params;
    fits.var2f.stats = stats;
    
    objresults(graph) = fits;
    save 'objective_fits' objresults results graphs;
end

