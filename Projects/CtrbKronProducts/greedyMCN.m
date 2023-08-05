function [C, ctrls] = greedyMCN(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = size(A,1);
ctrbs = cell(n,1);
for i=1:n
    Bi = zeros(n,1);
    Bi(i) = 1;
    ctrbs{i} = ctrb(A,Bi);
end

C = [];
ctrls = [];
D = 1:n;
while rank(C) < n
    for i=1:n
        if D(i) == 0
            continue;
        end
        maxi = 0;
        ri = rank([C ctrbs{i}]);
        if ri > maxi
            maxi = ri;
            argi = i;
        end
    end
    D(argi) = 0;
    C = [C ctrbs{argi}];
    ctrls = [ctrls argi];
end

end

