function [O] = getObservabilityMatricesNumeric(HG, x)
%GETOBSERVABILITYMATRICES This function computes the observability matrices
% for every vertex given a specified state x.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 2023

A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
n = size(A, 1);             % Get n
k = length(size(A));

% Unfold A
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix

% Calculate Middle part of expressions C (A B2 ... Bn) x
M = cell(n,1);
M{1} = eye(n);
M{2} = Amat;
for i=3:n
    M{i} = M{i-1} * getBp(sparse(Amat), i-1, k);
end

% Calculate numerical tradient
% Comment out: x = rand(n,1)
Mx = cell(n,1);
for i=1:n
    Mx{i} = M{i} * vecPower2(x, size(M{i}, 2));
end

delta = 1e-8;
gradM = cell(1,n); % Right part of expressions grad C (A B2 ... Bn x) for all xi
for j=1:n % Term in [eye(n), A, AB2, AB2B3, ..., AB2...Bn]
    % Ri = cell(n,1); % Right part of expressions C (A B2 ... Bn x) for specific i
    gradMj = zeros(n,n);
    for vx=1:n
        xh = x; xh(vx) = xh(vx) + delta;
        xl = x; xl(vx) = xl(vx) - delta;

        % Compute upper and lower steps
        h = M{j} * vecPower2(xh, size(M{j}, 2));
        l = M{j} * vecPower2(xl, size(M{j}, 2));
        
        % Compute the gradient of term M(j) with respect to vertex vx
        gradMj(vx,:) = (h - l) / (2 * eps);
    end
    gradM{j} = gradMj;
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

