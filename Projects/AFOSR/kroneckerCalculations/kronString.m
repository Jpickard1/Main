function S = kronString(S1, S2)
%KRONSTIRNG
%
%   Computes the kronecker product between 2 string arrays.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 27, 2023

S = string(zeros(length(S1) * length(S2), 1));
for i=1:length(S1)
    for j=1:length(S2)
        S((i-1)*length(S2) + j) = "(" + S1(i) + "" + S2(j) + ")";
    end
end

end