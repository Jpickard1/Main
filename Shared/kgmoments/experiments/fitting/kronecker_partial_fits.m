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
graphs = {'ca-GrQc','ca-HepPh','ca-HepTh','as20000102'};

%% Iterative over graphs
% For each graph, we want to compute a fit using both approaches: direct
% minimization and grid search.  

fitresults = containers.Map;

for gi=1:length(graphs)
    graph = graphs{gi};
    gdata = results(graph);
    
    fits = [];
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50);
    
    fits.all.params = params;
    fits.all.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'notris',true);
    
    fits.notris.params = params;
    fits.notris.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'notripins',true);
    
    fits.notripins.params = params;
    fits.notripins.stats = stats;
    
    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'nowedges',true);
    
    fits.nowedges.params = params;
    fits.nowedges.stats = stats;

    [params,stats] = kron_moment_fit(gdata,'alg','all','nstarts',50,'noedges',true);
    
    fits.noedges.params = params;
    fits.noedges.stats = stats;
    
    fitresults(graph) = fits;
    save 'partial_feature_fits' fitresults results graphs;
end

