function [xmin,fmin]=gridmin(f, X, verbose)

fmin=Inf;
xmin=NaN*ones(size(X,1),1);
if ~exist('verbose','var') || isempty(verbose), verbose=1000; end

for i=1:size(X,2)
    pts = X(:,i);
    
    fval = f(pts);
    if fval<fmin, fmin=fval; xmin=pts; end
    
    if mod(i,verbose)==0, fprintf('%7i : fmin=%20g\n', i, fmin); end
end