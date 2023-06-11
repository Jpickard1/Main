function A = rmat2_impl(P,r,nedges,noise)
% RMAT : generate power-law directed graph with R-MAT algorithm
%
% A = rmat(P,r)   returns a graph (adj matrix) with 2^r vertices,
%                 with vertex numbers not randomized.
% Implementation of the Recursive MATrix (R-MAT) power-law graph
% generation algorithm (Chakrabati, Zhan & Faloutsos).

% Original by Jeremy Kepner, repackaged by John R. Gilbert, summer 2006
% Modified by David Gleich, 2009-02-03, and again on 2009-02-13 to
% implement the algorithm from the paper

% Set number of vertices.
scale = r;
lgNv = scale;
Nv = size(P,1)^lgNv;

% % Set R-MAT probabilities.
% % Create a single parameter family.
% % Parameters can't be symmetric in order
% % for it to be a power law distribution.
% if ~exist('nedges','var') || isempty(nedges)
%     nedges = ceil(((sum(sum(P)))^scale - (sum(diag(P)))^scale));
% end

% P = P./sum(sum(P));

if nedges==0, A = logical(sparse(Nv,Nv)); return; end
P = P + noise*rand(size(P,1));
pprob = P(:)./sum(sum(P));
ncells = length(pprob);

if r==0
    A = sparse(nedges>=1);
else
    vals=rand(nedges,1);
    RE=zeros(size(P));
    for i=1:nedges
        for j=1:ncells
            if vals(i)<=pprob(j),
                RE(j)=RE(j)+1;
                continue
            end
        end
    end
    if r==1 % skip one recursive step
        A=sparse(RE>0);
    else
        A=cell2mat(arrayfun(@(x) rmat2_impl(P,r-1,x,noise),...
                RE, 'UniformOutput', false));
    end    
end
    