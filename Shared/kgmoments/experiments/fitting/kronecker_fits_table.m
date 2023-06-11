function kronecker_fits_table
%% Display a tabular of the fitting results from the Kronecker fits
% Written as a function so we can use inline functions

% The headers of the table 
%Graph	Type	Parameters			Features                                        Time
%               a	b	c           Vertices	Edges	Wedges	Tripins	Triangles

leskovec_kronfit_results = leskovec_fits;
load graphs_fits



for gi=1:length(graphs)
    graph = graphs{gi};
    fits = fitresults(graph);
    gdata = results(graph);
    
    r = ceil(log2(gdata.nverts));
    
    gname = graph;
    
    % output one line of data for the source
    output_row(gname,'Source', ...
        '---', ... %a, b, c
        gdata, '---', gdata.compute_time);
    
        gname = ''; % skip subsequent output of gname
    
    output_row(gname,'Direct', ...
        [fits.direct.params.a, fits.direct.params.b, fits.direct.params.c], ...
        expected_kronecker_moments(fits.direct.params), ...
        fits.direct.stats.obj, fits.direct.stats.fit_time, gdata,0);

    output_row(gname,'Grid', ...
        [fits.grid.params.a, fits.grid.params.b, fits.grid.params.c], ...
        expected_kronecker_moments(fits.grid.params), ...
        fits.grid.stats.obj, fits.grid.stats.fit_time, gdata,1);    
    
    if all(fits.first.params.mat~=0)
        output_row(gname,'Leading', ...
            [fits.first.params.a, fits.first.params.b, fits.first.params.c], ...
            expected_kronecker_moments(fits.first.params), ...
            fits.first.stats.obj, fits.first.stats.fit_time, gdata,1);
    end
       
    
    if isKey(leskovec_kronfit_results,graph)
        abc = leskovec_kronfit_results(graph);
        
        % copied from kron_moment_fit
        weight = [1,1,1,1];
        f = [];
        f(1) = gdata.nedges;
        f(2) = gdata.nwedges;
        f(3) = gdata.ntripins;
        f(4) = gdata.ntris;
        f = f';
        
        e = @(a,b,c,r) [
            (1/2)*((a+2*b+c)^r - (a+c)^r)
            (1/2)*(((a+b)^2 + (b+c)^2)^r - 2*(a*(a+b)+c*(c+b))^r - (a^2 + 2*b^2+c^2)^r + 2*(a^2 + c^2)^r)
            (1/6)*(((a+b)^3 + (b+c)^3)^r - 3*(a*(a+b)^2 + c*(b+c)^2)^r - 3*(a^3 + c^3 + b*(a^2+c^2)+b^2*(a+c) + 2*b^3)^r + 2*(a^3 + 2*b^3 + c^3)^r + 5*(a^3 + c^3 + b^2*(a+c))^r + 4*(a^3 + c^3 + b*(a^2 + c^2))^r - 6*(a^3 + c^3)^r)
            (1/6)*((a^3 + 3*b^2*(a+c) + c^3)^r - 3*(a*(a^2 + b^2) + c*(b^2 + c^2))^r + 2*(a^3+c^3)^r)];

        curobj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r))./(max(f,eps(1)) ));
        
        output_row(gname,'KronFit',...
            abc, ...
            expected_kronecker_moments(abc,r), curobj(abc), '---', gdata,1);
        % get KronFits!
    end
    
    fprintf('\\midrule \n');
end

