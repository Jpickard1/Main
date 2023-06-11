function [rmat_params,stats,G] = kron_moment_fit(gdata,varargin)
% KRON_MOMENT_FIT Fit the parameter of a 2x2 Kronecker model by matching moments
%
% The moments of an adjacency matrix can be used to count features of a
% graph, such as the number of edges, number of wedges, number of tripins,
% and number of triangles.  There are efficient algorithms for each of
% these moment computations.  Matching them is a simple and scalable way to
% find reasonable parameters for a Kronecker model.
%
% params = kron_moment_fit(A) takes an adjacency matrix A and computes the
% moments of the adjacency matrix using the MatlabBGL library via the
% function GRAPH_PROPERTIES.  These are then used in fitting procedure.
% Alternatively,
% ... = kron_moment_fit(gdata) 
% takes a structure of data about the graph where:
%     gdata.nverts: the number of vertices
%     gdata.nedges: the number of undirected edges (nnz(A)/2)
%     gdata.nwedges: the number of wedges
%     gdata.ntripins: the number of tripins
%     gdata.ntris: the number of triangles (optional, see below)
% OR gdata can be the degree distribution for the graph.  In the latter
% case, the fit is computed by matching the moments for the number of
% edges, the number of wedges, and number of tripins, which can all be
% computed just from the degree distribution.
%
% [params,stats] = kron_moment_fit(gdata) returns additional data about the
% fitting procedure, including the time for the fitting algorithm, and the
% final objective function.
%
% [params,stats,G] = kron_moment_fit(gdata) also produces a realization of
% the graph from the parameters.  This version of the command outputs a set
% of statistics on the fit, unless this is overridden below.
%
% ... = kron_moment_fit(gdata,'param1',value1,'param2',value2,...) changes 
% a few optional configuration settings.  The possible parameters are:
%
%  'obj' -- [{'var2'} | 'var' | 'diff1']
%    a string parameter that changes the objective function for the
%    parameter choice.  The objectives are:
%      'var2' = sum_F (F-E(F|a,b,c))^2/E(F|a,b,c)^2
%      'var' = sum_F (F-E(F|a,b,c))^2/E(F|a,b,c)
%      'diff1' = sum_F |F-E(F|a,b,c)|/|F|
%    Using 'var' is best when you expect a good fit.  Using 'diff1' is best
%    when you do not expect a good fit.  Using var2 is a smooth version of
%    'diff1' and is a good general choice.
%    
%  'alg' -- [{'direct'} | 'grid' | 'first' | 'all']
%    a string parameter indicating the algorithm to use to fit the 
%    moments of the distribution.  The direct method uses Matlab's 
%    fmincon.  The grid method divides the search region into a grid
%    and does exhaustive search.  The 'first' method matches the first
%    order terms in the moment expressions.  The 'all' approach is to run
%    all algorithms and report the minimum.
%
%  'grid_size' -- [{100} | positive integer]
%    a numeric parameter that gives the grid size to use in the search.
%    The total number of grid points is grid_size^3 because there are three
%    parameters.  The grid search is written using parfor, and thus can
%    take advantage of multiple Matlab's via this coarse parallelization.
%
%  'nstarts' -- [{10} | positive integer]
%    a numeric parameter for the number of random starting points to use
%    for the direct search method.
%
%  'sampletype' -- [{'rmat'} | 'coinflip' | 'skron' | 'rmatn']
%    a string parameter naming the procedure to use to generate a sample of
%    the fitted parameters if the output G is requested.  Both 'coinflip'
%    and 'skron' call the skrongraph function.  The rmat option uses the
%    rmat function and the rmatn option uses the rmatn function with the
%    maximum possible noise.
%
%  'notris' -- [{false} | true]
%    a logical parameter indicating that the triangle count is
%    not included in the features to fit.  This is automatically 
%    enabled if gdata does not contain triangles or if only a 
%    degree vector is passed.
%
%  'notripins' -- [{false} | true]
%    a logical parameter indicating that the tripin count is
%    not included in the features to fit.  
%
%  'nowedges' -- [{false} | true]
%    a logical parameter indicating that the wedge count is
%    not included in the features to fit.  
%
%  'noedges' -- [[{false} | true]
%    a logical parameter indicating that the edge count is
%    not included in the features to fit.  
%
%  'display' -- [{true} | false]
%    a logical parameter that determines whether or not properties of the
%    fit are reported by the algorithm.  If display is false, then the
%    properties of the sampled graph G are _not_ computed.
%
%  'sampleprops' -- [{false} | true]
%    generate properties from a sample of the parameters even if the sample
%    output is not requested.
%  
% Example:
%   A = rmat([1,0.7;0.7,0.3],10);
%   props = kron_moment_fit(A)


optsu = struct(varargin{:});
opts = struct('obj','var2f','alg','direct','grid_size',100,'sampletype','rmat',...
    'noedges', false, 'notris',false,'notripins',false,'nowedges',false,...
    'display',true, 'nstarts', 10, 'sampleprops', false);
for f=fieldnames(optsu)', fs=f{1}; opts.(fs) = optsu.(fs); end

if isstruct(gdata)
    if ~isfield(gdata,'ntris')
        if opts.display
            fprintf('Disabling triangle moments because ntris not in gdata\n');
        end
        opts.notris = true;
    end
else
    if size(gdata,1) == size(gdata,2)
        % this must be an adjacency matrix
        gdata = graph_properties(gdata);
    else
        % this is a degree distribution
        % TODO abstract this distribution
        t0 = tic;
        gdata.nverts = length(degs);
        gdata.nedges = sum(degs)/2;
        gdata.nwedges = sum(d.*(d-1))/2;
        gdata.compute_time = toc(t0);
    end
end

stats.gdata = gdata;
stats.moment_compute_time = gdata.compute_time;
    
% Setup expectations
r = ceil(log2(gdata.nverts));
e = @(a,b,c,r) [
    (1/2)*((a+2*b+c)^r - (a+c)^r)
    (1/2)*(((a+b)^2 + (b+c)^2)^r - 2*(a*(a+b)+c*(c+b))^r - (a^2 + 2*b^2+c^2)^r + 2*(a^2 + c^2)^r)
    (1/6)*(((a+b)^3 + (b+c)^3)^r - 3*(a*(a+b)^2 + c*(b+c)^2)^r - 3*(a^3 + c^3 + b*(a^2+c^2)+b^2*(a+c) + 2*b^3)^r + 2*(a^3 + 2*b^3 + c^3)^r + 5*(a^3 + c^3 + b^2*(a+c))^r + 4*(a^3 + c^3 + b*(a^2 + c^2))^r - 6*(a^3 + c^3)^r)
    (1/6)*((a^3 + 3*b^2*(a+c) + c^3)^r - 3*(a*(a^2 + b^2) + c*(b^2 + c^2))^r + 2*(a^3+c^3)^r)];

f = zeros(4,1);
f(1) = gdata.nedges;
f(2) = gdata.nwedges;
f(3) = gdata.ntripins;
f(4) = gdata.ntris;

%% Minimize
t0 = tic;
if opts.notris
    weight=[1,1,1,0];
elseif opts.notripins    
    weight=[1,1,0,1];
elseif opts.nowedges
    weight=[1,0,1,1];
elseif opts.noedges
    weight=[0,1,1,1];
else
    weight=[1,1,1,1];
end

switch opts.obj
    case {'var2','var2e'}
        obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r)).^2./(max(e(x(1),x(2),x(3),r).^2,eps(1)) ));
    case {'var','vare'}
        obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r)).^2./(max(e(x(1),x(2),x(3),r),eps(1)) ));
    case {'varf'}
        obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r)).^2./(max(f,eps(1)) ));
    case {'var2f'}
        obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r)).^2./(max(f.^2,eps(1)) ));
    case {'diff1','diff1f'}
        obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r))./(max(f,eps(1)) ));
    case {'diff1var','diff1e'}
        obj = @(x) sum( diag(weight)*abs(f - e(x(1),x(2),x(3),r))./(max(e(x(1),x(2),x(3),r),eps(1)) ));
    otherwise
        error('kron_moment_fit:invalidParameter',...
            'the objective ''%s'' is not recognized', opts.obj);
end

minobj = Inf;
xmin = [];

if strcmp(opts.alg,'direct') || strcmp(opts.alg,'all')
    confun = @(x) deal( x(3)-x(1) , [] ); % c-a <= 0 => c<=a
    lb = [0 0 0];
    ub = [1 1 1];
    options = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Algorithm','active-set','Display','off');

    for i=1:opts.nstarts
        x0 = sort(rand(1,3),2,'descend');
        x = fmincon(obj,x0,[],[],[],[],lb,ub,confun,options);
        % enforce constraints exactly!
        x = max(x,0);
        x = min(x,1);
        %[x0; x]
        if obj(x)<minobj, xmin = x; minobj = obj(x); end
    end
    stats.direct.params = x;
    stats.direct.obj = obj(x);
end

if strcmp(opts.alg,'grid') || strcmp(opts.alg,'all')
    npts = opts.grid_size;
    [x1,x2,x3]=ndgrid(0:1/npts:1, 0:1/npts:1, 0:1/npts:1);
    X=[x1(:), x2(:), x3(:)]';
    X(:,x1(:) < x3(:)) = [];
    x=gridmin(obj, X, 0);
    stats.grid.params = x;
    stats.grid.obj = obj(x);
    if obj(x)<minobj, xmin = x; minobj = obj(x); end
end

if ~opts.notris && ~opts.nowedges && ~opts.notripins && ...
        (strcmp(opts.alg,'first') || strcmp(opts.alg,'all'))
    first_e = (2*gdata.nedges)^(1/r); % using notation from the paper
    first_h = (2*gdata.nwedges)^(1/r);
    first_tris = (6*gdata.ntris)^(1/r); % the tripins are redundant in this one
    
    if 2*first_h - first_e^2 < 0 
        if strcmp(opts.alg,'first')
            warning('kron_moment_fit:unknownFit',...
                ['matching leading terms is only possible when\n' ...
                'var(degs) >= mean(degs)']);
            x = [0 0 0];
        else
            % in this case, they are doing other algorithms too,
            % so we don't need to worry
        end
    else
        x = 0.5*(first_e + sqrt(2*first_h - first_e^2));
        y = 0.5*(first_e - sqrt(2*first_h - first_e^2));
        % these should be true for any graph
        assert(x >= y);
        assert(y >= 0);
        bmin = max(x-1,y-1);
        bmin = max(bmin, 0);
        bmax = min(x,y); % 
        bmax = min(bmax,1);
        % make sure the region isn't empty
        if bmin <= bmax
            uniobj = @(b) abs( ...
                (x-b).^3 + (y-b).^3 + 3.*b.^2.*((x-b) + (y-b)) ... % expected triangle count
                - first_tris);
            bs = linspace(bmin,bmax,10000);
            objs = uniobj(bs);
            
            [~,index] = min(objs);
            b = bs(index);
            x = [x-b b y-b];
        else
            warning('kron_moment_fit:unknownFit',...
                ['matching leading is not possible, using best approximation']);
            x = [0 0 0];
%             uniobj = @(b) abs( ...
%                 (min(x-b,1)).^3 + (max(y-b,0)).^3 + ...
%                     3.*b.^2.*(min(x-b,1) + max(y-b,0)) ... % expected triangle count
%                 - first_tris);
%             bs = linspace(0,1,10000);
%             objs = uniobj(bs);
%             if min(objs) < minobj
%                 [~,index] = min(objs);
%                 b = bs(index);
%                 xmin = [min(x-b,1) b max(y-b,0)];
%                 assert(obj(xmin) < minobj);
%                 minobj = obj(xmin);
%             end
        end
    end
    stats.leading.params = x;
    stats.leading.obj = obj(x);
    if obj(x) < minobj
        xmin = x;
        minobj = obj(xmin);
    end
end



stats.fit_time = toc(t0);

%% output
x = xmin;
a = x(1);
b = x(2);
c = x(3);
rmat_params=struct();
rmat_params.a = a;
rmat_params.b = b;
rmat_params.c = c;
rmat_params.r = r;
rmat_params.mat = [a b; b c];

stats.obj = obj(x);

fitdata = e(a,b,c,r);
stats.fitdata.nedges = fitdata(1);
stats.fitdata.nwedges = fitdata(2);
stats.fitdata.ntripins = fitdata(3);
stats.fitdata.ntris = fitdata(4);

%% report
if opts.display
    fprintf('\n');
    fprintf('Kronecker Moment Fit\n');
    if isfield(gdata,'name')
        fprintf('graph               = %s\n', gdata.name);
    end
    fprintf('src vertices        = %20g\n', gdata.nverts);
    fprintf('objective value     = %20g\n', obj(x));
    fprintf('fitted params a,b,c = [%6.4f , %6.4f , %6.4f]\n', a,b,c);
    
    fprintf('%5s',''); 
     fprintf(' %17s %17s %17s %17s','edges','hairpins','tripins','triangles');
     fprintf('\n');
    fprintf('%5s','src'); fprintf(' %17i', f); fprintf('\n');
    fprintf('%5s','est'); fprintf(' %17i', round(e(a,b,c,r))); fprintf('\n');
    
    if nargout>2 || opts.sampleprops
        t0 = tic;
        switch opts.sampletype
            case 'rmat'
                G = rmat([a b; b c], r);
            case 'rmatn'
                G = rmatn([a b; b c], r, 1.);
            case {'coinflip','skron'}
                G = skrongraph([a b; b c], r);
            otherwise
                error('kron_moment_fit:invalidParameter',...
                    'unknown sampletype ''%s''', opts.sampletype);
        end
        stats.sample_compute_time = toc(t0);
            
        G = triu(G,1); G=G|G';
        G = G - diag(diag(G));
        props = graph_properties(G);
        g(1) = props.nedges;
        g(2) = props.nwedges;
        g(3) = props.ntripins;
        g(4) = props.ntris;
        stats.sample_gdata = props;
        fprintf('%5s','fit'); fprintf(' %17i', g); fprintf('\n');
    end
end

