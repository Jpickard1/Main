function output_row(graph,type,abc,gdata,obj,dt,relative,dupverts)

% Looking closely at the tables, I think I see
% a better way to present the information.  Hopefully
% it is just a short programming change.
% 
% Tables 2:4
% 
% objective in scientific notation
% source counts as integers
% expected counts .. as a ratio to the source count
% 
% some of the ratios will have to be in
% scientific notation.  In other cases if
% they're all close to 1 then maybe a fixed
% point representation would be better.
% I think we need to be consistent within
% blocks defined by a specific (graph, feature) pair.
% For non-scientific notation, it's best to have
% the decimal points line up.  There is a fussy
% way to do that in LaTeX.  I usually get by
% via right justifying those columns and if
% necessary putting in some ad hoc space.

% dupverts = 1 if we should output a duplicate col token
if ~exist('dupverts','var') || isempty(dupverts),
    dupverts = 0;
end


if ~isempty(graph)
    fprintf('\\multicolumn{2}{l}{ \\bfseries %20s} \\\\ \n', graph);
end
    
fprintf('%9s ', type);

if ischar(abc)
    fprintf('& %6s & %6s & %6s ', '---','---','---'); % no a,b,c for source
else
    fprintf('& %6.3f & %6.3f & %6.3f ', abc(1), abc(2), abc(3));
end

if exist('relative','var') && ~isempty(relative)
    
    if ~dupverts
        fprintf('& %15ld ', gdata.nverts);
    else
        fprintf('& %15s ', '\dupcolval');
    end
    
    fprintf('& %15.2f ', gdata.nedges/relative.nedges)
    fprintf('& %15.2f ', gdata.nwedges/relative.nwedges)
    fprintf('& %15.3f ', gdata.ntripins/relative.ntripins)
    fprintf('& %15.4f ', gdata.ntris/relative.ntris)
   
        
else
    nums = round([gdata.nverts, gdata.nedges, gdata.nwedges, gdata.ntripins, gdata.ntris]);
    strs = format_numbers(nums);
    for i=1:length(strs)
        fprintf('& %15s ', strs{i});
    end
end


if ischar(obj)
    fprintf('& %3s ', obj);
else
    %objs = latex_sci(obj);
    %fprintf('& %15s ', objs{1});
    fprintf('& %5.3f ', obj);
end

if ~isempty(dt)
    if ischar(dt)
        fprintf('& %6s ', dt);
    else
        if dt>0.05
            fprintf('& %6.1f', dt);
        else
            fprintf('& $<\\!\\!0.05$');
        end
    end
end

fprintf('\\\\ \n');

function strs=format_numbers(nums)
% Make the number nicely formatted for latex

strs = {};
for n=nums(:)'
    if n<1e8 && n>-1e8
        numtxt=num2str(n,'%9i');
    else
        numtxt = latex_sci(n);
        numtxt = numtxt{1};
    end
    strs{end+1} = numtxt;
end

