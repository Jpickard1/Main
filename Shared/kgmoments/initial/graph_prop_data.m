addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl/4.0');

results = containers.Map;

%%

graphs = {
  '/home/dgleich/data/graph-db/dgleich/wb-cs.stanford.smat'
  '/home/dgleich/data/graph-db/misc/itdk0304.smat'
  '/home/dgleich/data/graph-db/dgleich/usroads-cc.smat'
  '/home/dgleich/data/graph-db/snap/ca-HepTh.smat'
  '/home/dgleich/data/graph-db/snap/cit-HepTh.smat'
  '/home/dgleich/data/graph-db/snap/p2p-Gnutella31.smat'
  '/home/dgleich/data/graph-db/dgleich/wikipedia-20051105.smat'
};

%%

for gi=1:length(graphs)
    g = graphs{gi};
    [path,name,ext]=fileparts(g);
    name
    if 
    A = readSMAT(g);
    props = graph_properties(A,name);
    
    results(name) = props;
end
save 'graph_props.mat' results;
    