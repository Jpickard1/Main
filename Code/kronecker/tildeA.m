function [A] = tildeA(A, mb, nb)
%TILDEA Set up the tilde matrix T from A subject to the dimension of B
% The dimensions of A must be divisible by the corresponding dimensions of
% B
%
%   This code was taken directly from the thesis of Nikos P. Pitasianis
%   titled: The Kronecker Product In Approximation and Fast Transform
%   Generation
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 3, 2023

[m, n] = size(A);
mc = m / mb; nc = n / nb;
T = zeros(mb * nb, mc * nc);
x = zeros(1, mc * nc);      % as a block stretcher

for ib=1:mb
    for jb=1:nb
        x(:) = A((ib-1)*mc+1:ib*mc,(jb-1)*nc+1:jb*nc); % stretch bock
        T((jp-1)*mb+ib,:) = x;
    end
end

end

