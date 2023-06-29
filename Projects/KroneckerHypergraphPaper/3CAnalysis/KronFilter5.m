function [B, C] = KronFilter5(A, b, c)
%KRONFILTER5 A has blocks of size b x c
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 28, 2023

[m, n] = size(A);

B = zeros(round(m/b), round(n/c));
C = zeros(b,c);

for i=1:(m/b)
    for j=1:(n/c)
        B(i,j) = sum(A((i-1)*b+1:i*b, (j-1)*c+1:j*c), 'all');
        if B(i,j) ~= 0
            C = C + (A((i-1)*b+1:i*b, (j-1)*c+1:j*c) / (B(i,j)));
        else
            disp('skip')
        end
    end
end
C = C ./ (i*j);

end

