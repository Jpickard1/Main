function objectives_table
%% Display a tabular of the fitting results from the Kronecker fits
% Written as a function so we can use inline functions

% The headers of the table 
%Graph	Type	Parameters			Features                                        Time
%               a	b	c           Vertices	Edges	Wedges	Tripins	Triangles

load objective_fits

for gi=1:length(graphs)
    graph = graphs{gi};
    fits = objresults(graph);
    gdata = results(graph);
    
    r = ceil(log2(gdata.nverts));
    
    gname = graph;
    
    % output one line of data for the source
    output_row(gname,'Source', ...
        '---', ... %a, b, c
        gdata, '---',[],[]);
    
    gname = ''; % skip subsequent output of gname
        
    handle_row(gname,'$D_{\text{sq}}, N_{\e}$',fits.var, gdata, 0);
    handle_row(gname,'$D_{\text{sq}}, N_{\e^2}$',fits.var2, gdata);
    handle_row(gname,'$D_{\text{sq}}, N_{F}$',fits.varf, gdata);
    handle_row(gname,'$D_{\text{sq}}, N_{F^2}$',fits.var2f, gdata);
    handle_row(gname,'$D_{\text{abs}}, N_{\e}$',fits.diff1e, gdata);
    handle_row(gname,'$D_{\text{abs}}, N_{F}$',fits.diff1, gdata);
    
    % see if any of them get a leading fit (they all would have...)
    if isfield(fits.var.stats,'leading')
        params = fits.var.stats.leading.params;
        fits.var.stats.leading.params = struct('a',params(1),'b',params(2),'c',params(3),'r',r);
        fits.var.stats.leading.stats.obj = '---';
        handle_row(gname,'Leading',fits.var.stats.leading,gdata);
    end
       
    fprintf('\\midrule \n');
end

function handle_row(gname,label,fit,gdata, dup)
   
if ~exist('dup','var'), dup=1; end

if ~ischar(fit.stats.obj)
    objstr = latex_sci(fit.stats.obj);
    fit.stats.obj = objstr{1};
end
    
     output_row(gname,label, ...
        [fit.params.a, fit.params.b, fit.params.c], ...
        expected_kronecker_moments(fit.params), ...
        fit.stats.obj, [], gdata, dup);

