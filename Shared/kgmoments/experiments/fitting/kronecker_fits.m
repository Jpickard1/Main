%% Fit Kronecker Parameters to real-world networks via moment matching
% This experiment tests the kron_moment_fit routine to determine the
% parameters that it produces for real-world networks.

%% Setup

addpath('../../matlab');
addpath('~/dev/matlab-bgl');

%% Load precomputed graph feature data
% This experiment uses precomputed graph feature data exclusively.
% To add graphs to the experiment, add graphs to the
% 'data/graph_features.m' script.  

load('../../data/graphs_props.mat');

%%
% The output of loading this script is a containers.Map type called results
graphs = results.keys

%% Iterative over graphs
% For each graph, we want to compute a fit using both approaches: direct
% minimization and grid search.  

fitresults = containers.Map;

for gi=1:length(graphs)
    graph = graphs{gi};
    gdata = results(graph);
    
    fits = [];
    
    [params,stats] = kron_moment_fit(gdata,'alg','direct','nstarts',50);
    
    fits.direct.params = params;
    fits.direct.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','grid');
    fits.grid.params = params;
    fits.grid.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','first');
    fits.first.params = params;
    fits.first.stats = stats;
    
    fitresults(graph) = fits;
    save 'graphs_fits' fitresults results graphs;
end

