function [N1, NE, NI, NP] = tensorNorms(T)
%TENSORNORMS Computes norm of a multi-dimensional array
%
%   N1 = l1 norm
%   NE = Euclidean norm
%   NI = infintie or max norm
%   NP = inner product norm
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 2, 2023

% 1 norm
S = sum(T);
while numel(S) > 1
    S = sum(S);
end
N1 = S;
% Euclidean norm
E = T.^2;
S = sum(E);
while numel(S) > 1
    S = sum(S);
end
NE = sqrt(S);
% Infity norm
S = max(T);
while numel(S) > 1
    S = max(S);
end
NI = S;

% Inner product with itself
I = T.^2;
S = sum(I);
while numel(S) > 1
    S = sum(S);
end
NP = S;

end