function [ll,grad] = kronloglike(A,T,warmup,samples)
% KRONLOGLIKE Estimate the log-liklihood of a Kronecker fit using SNAP
%
% [ll,grad] = kronloglike(A,T)
%
% TODO -- document

[i,j] = find(A);
[mydir myfilename] = fileparts(mfilename('fullpath'));
loglikeprog = fullfile(mydir,'kronloglike');
graphfile = fullfile(mydir,'matlab-graph.txt');
f = fopen(graphfile,'w');
for k=1:nnz(A)
    fprintf(f,'%i %i\n',i(k),j(k));
end
fclose(f);

sizearg = sprintf('-n:%i',size(T,1));
%mattext = sprintf(' %f',T(:));
mattext = mat2str(T);
mattext = mattext(2:end-1);% remove []'s
matarg = sprintf('-m:"%s"',mattext);
grapharg = sprintf('-i:%s',graphfile);
warmuparg = sprintf('-w:%i',round(warmup));
samplesarg = sprintf('-s:%i',round(samples));

args = sprintf('%s %s %s %s %s',grapharg,sizearg,matarg,warmuparg,samplesarg);

cmd=sprintf('%s %s',loglikeprog,args);
cmd
system(cmd);
