function [O] = HGObsvNum(HG, x)
%HGOBSVSYM This funciton computes numeric observability matrices for
%   hypergraphs.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 20, 2023

delta = 1e-8;
n = size(HG.IM,1);

% Compute Jp vectors
gradM = cell(1,n); % Right part of expressions grad C (A B2 ... Bn x) for all xi
for i=n:-1:1
    gradMj = zeros(n,n);
    for vx=1:n
        xh = x; xh(vx) = xh(vx) + delta;
        xl = x; xl(vx) = xl(vx) - delta;
        h = Jp(HG, i, xh);
        l = Jp(HG, i, xl);
        gradMj(vx,:) = (h - l) / (2 * eps);
    end
    gradM{i} = gradMj;
    disp(i)
end

O = cell(n,1);
for vx=1:n
    Ci = zeros(1,n); Ci(vx) = 1;% Equation 9
    Ovx = zeros(n,n);
    for i=1:n
        Ovx(i,:) = Ci * gradM{i};
    end
    O{vx} = Ovx;
end

end