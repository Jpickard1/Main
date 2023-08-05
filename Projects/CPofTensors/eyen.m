function [I] = eyen(n,k)
%EYEN Creates the k-mode n-dimensional identity tensor
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2023

I = zeros(n * ones(1,k));
idxs = "i";
for i=1:k-1
    idxs = idxs + ",i";
end
cmd = "I(" + idxs + ")=1";
for i=1:n
    eval(cmd);
end

end

