function [A] = skrongraph(T,r,varargin)
% SKRONGRAPH Generate an instance of a Stochastic Kronecker graph
%
% G = skrongraph(T,r) based on the flipping a coin for each of the cells in
% the matrix.  This procedure is rather expensive due to the large number
% of coin flips that generate a tails result

optsu = struct(varargin{:});
opts = struct('symmetric',0,'noise',0);
for f=fieldnames(optsu)', f=f{1}; opts.(f) = optsu.(f); end

% check parameters
if ~isequal(size(T,1),size(T,2))
    error('The generator matrix must be square');
end

% the size of the initatior matrix
n = size(T,1);

% the size of the overall Kronecker graph
N = n^r;

% keep things in log-space for precision
lT = log(T);
% size of the tensor structure
siz = n*ones(1,r);
rowpos = zeros(r,1);
colpos = zeros(r,1);
ind2sub_k = [1 cumprod(siz(1:end-1))];


% keep a buffer of random data
bufsize = min(1e6,N^2);
rbuf = log(rand(bufsize,1));
bufpos = 1;

% keep a buffer of edges, over-allocated over the expectation, a little
edgebuf_size = ceil(1.1*((sum(sum(T)))^r - (sum(diag(T)))^r));
edgebuf = zeros(edgebuf_size,2);
edgepos = 1;



% for testing...
%P = zeros(N,N);
%A = zeros(N,N);


if opts.symmetric
    error('not implemented');
else    
    for i=1:N
        for j=1:N
            % determine the coordinates of the cell
            
            % call ind2sub here (in-place without a function call)
            % [rowpos{:}] = ind2sub(siz,i);
            ndx = i;
            for ind = r:-1:1,
                vi = rem(ndx-1, ind2sub_k(ind)) + 1;         
                vj = (ndx - vi)/ind2sub_k(ind) + 1; 
                rowpos(ind) = vj; 
                ndx = vi;     
            end
            
            % [colpos{:}] = ind2sub(siz,j);
            ndx = j;
            for ind = r:-1:1,
                vi = rem(ndx-1, ind2sub_k(ind)) + 1;         
                vj = (ndx - vi)/ind2sub_k(ind) + 1; 
                colpos(ind) = vj; 
                ndx = vi;     
            end
            
            logprob = 0;
            for k=1:r
                logprob = logprob + lT(rowpos(k),colpos(k));
            end
            % check if we'll get a cell here
            if bufpos > bufsize
                rbuf = log(rand(bufsize,1))+rand(1,1);
                bufpos = 1;
            end
            randval = rbuf(bufpos);
            bufpos = bufpos + 1;
            
            %P(i,j) = logprob;
            if randval < logprob
                if edgepos>edgebuf_size
                    edgebuf = [edgebuf; edgebuf]; %#ok<AGROW> % double the size
                    edgebuf_size = edgebuf_size*2;
                end
                edgebuf(edgepos,:) = [i,j];
                edgepos = edgepos + 1;
            end
        end
    end
end
edgebuf = edgebuf(1:edgepos-1,:); % resize to used area!
A = sparse(edgebuf(:,1),edgebuf(:,2),1,N,N);