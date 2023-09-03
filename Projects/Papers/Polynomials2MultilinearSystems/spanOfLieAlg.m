function [S] = spanOfLieAlg(p, B)
%SPANOFLIEALG Evaluates the span of the Lie algebra of a polynomial and the
% input signals near the origin.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 1, 2023

S = B;
i = 1;
while i <= size(S,2)
    b = S(:,i);
    fb = p.eval(b);

    % Check if the new term is in the image of S
    x = S \ fb;
    % if it is not in the image of S, append it to the image of S
    if sum(fb - S * x) ~= 0
        S = [S fb];
    end
    i = i + 1;
    if rank(S) == size(S)
        break;
    end
end

end

