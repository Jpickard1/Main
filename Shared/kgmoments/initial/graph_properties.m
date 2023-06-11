function props=graph_properties(A,name)

% normalize graph
A = A|A';
%A = largest_component(A);
A = A-diag(diag(A));
    
d = full(sum(A,2));

props = struct;
if nargin>1,
    props.name = name;
end
props.degs = d;
props.nverts = size(A,1);
props.nedges = sum(d)/2;
props.nwedges = sum(d.*(d-1))/2;
props.ntripins = sum(d.*(d-1).*(d-2))/6;
ccfs = clustering_coefficients(A,'undirected',1','unweighted',1);
tris = ccfs.*d.*(d-1);
props.ntris = sum(tris)/6;
props.clustcoeffs = ccfs;

cores = core_numbers(A);
props.cores = cores;
props.nonecore = sum(cores < 2);
