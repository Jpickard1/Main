function [A] = HG2Aten(HG)
%HG2ATEN Constructs Directed Adjacency Tensor
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 8, 2023

et1 = HG.E{1,1};
degree = length(et1{1});
n = numel(HG.V);
A = zeros(n * ones(1, degree + 1));
numE = size(HG.E,1);

for i=1:numE
    et = HG.E{i,1};
    eh = HG.E{i,2};
    if iscell(et)
        et = et{1};
    end
    if iscell(eh)
        eh = eh{1};
    end
    from = "";
    for j=1:length(et)
        from = from + "," + string(et(j));
    end
    for j=1:length(eh)
        cmd = "A(" + string(eh(j)) + from + ") = 1;";
        eval(cmd);
        disp(cmd);
    end
end

end

