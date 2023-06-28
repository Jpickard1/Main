function [B, C] = HyperKronFilter(A, Bs)
%HYPERKRONFILTER A fast approxixmation for A = B kron C where B and C are
% unknown 3-way tensors
%   
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 28, 2023

n = size(A,1);
B = zeros(floor(n / Bs), floor(n / Bs), floor(n / Bs));
C  = zeros(Bs, Bs, Bs);
for i=1:size(B,1)
    for j=1:size(B,1)
        for k=1:size(B,1)
            B(i,j,k) = mean(A((i-1)*Bs+1:i*Bs, (j-1)*Bs+1:j*Bs, (k-1)*Bs+1:k*Bs), 'all');
            if B(i,j,k) ~= 0
                C = C + (A((i-1)*Bs+1:i*Bs, (j-1)*Bs+1:j*Bs, (k-1)*Bs+1:k*Bs)); % / (B(i,j,k)/n^3));
            else
                disp('skip')
            end
        end
    end
end
% F = median(Fp,3);
% C = C ./ (numel(A));

end

