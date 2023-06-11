function props=graph_properties(A,name,notris)
% GRAPH_PROPERTIES Compute moments of the adjacency matrix

if ~exist('notris','var') || isempty(notris)
    notris = false;
end

% normalize graph
A = A|A';
A = A-diag(diag(A));
   
t0 = tic;
d = full(sum(A,2));

props = struct;
if nargin>1 && ~isempty(name)
    props.name = name;
end
props.degs = d;
props.nverts = size(A,1);
props.nedges = sum(d)/2;
props.nwedges = sum(d.*(d-1))/2;
props.ntripins = sum(d.*(d-1).*(d-2))/6;
if ~notris
    ccfs = clustering_coefficients(A,'undirected',1','unweighted',1);
    tris = ccfs.*d.*(d-1);
    props.ntris = sum(tris)/6;
    props.clustcoeffs = ccfs;
end
dt = toc(t0);
props.compute_time = dt;