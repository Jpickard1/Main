load 'graph_props.mat'

for name=results.keys
    name = name{1};
    [params,obj,G]= kron_moment_fit(results(name));
end

%%
% fitting without triangle counts
for name=results.keys
    name = name{1};
    [params,obj,G]= kron_moment_fit(results(name),1);
end