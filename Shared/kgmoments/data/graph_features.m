%% Graph features
%

addpath('~/dev/matlab');
addpath('~/dev/matlab-bgl');

%%

graphs = {
  'as20000102.smat'
  'ca-GrQc.smat'
  'ca-HepPh.smat'
  'ca-HepTh.smat'
  'as-22july06.smat'
  'usroads-cc.smat'
  'as-skitter.smat'
  'wikipedia-20051105.smat'
  'hollywood-2009.smat'
};
%%
results = containers.Map;

for gi=1:length(graphs)
    g = graphs{gi};
    [path,name,ext]=fileparts(g);
    name
    A = readSMAT(g);
    props = graph_properties(A,name);
    results(name) = props;
end
save 'graphs_props.mat' results;
    