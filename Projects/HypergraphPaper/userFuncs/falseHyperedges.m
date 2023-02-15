function [F] = falseHyperedges(HGtrue, numFalse)
%FALSEHYPEREDGES Generate false hyperedges for hyperedge prediction
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: January 29, 2023

% Params
alpha = 0.90;

[n, m] = size(HGtrue.IM);

F = cell(m, 1);
for i=1:numFalse
    % copy edge
    ei = randsample(1:m,1);
    % select ith hyperedes
    e = find(HGtrue.IM(:, ei) ~= 0);
    % if ei is a singleton continue
    if length(e) < 2
        continue
    end
    % Select number of vertices to replace
    numReplace = round(alpha * length(e));
    % Select which vertices will be replaced
    v = randsample(1:length(e), numReplace);
    % Get list of possible vertices to replace (i.e. vertices not in ei)
    replacements = setdiff(1:n, e);
    % Select each replacement
    for j=1:length(v)
        e(v(j)) = randsample(replacements, 1);
    end
    F{i} = e;
end

T = cell(m,1);
for i=1:m
    e = sort(find(HGtrue.IM(:,i) ~= 0));
    T{i} = char(join(string(e), "-"));
end

% Check false hyperedges are valid
%   1. not true hyperedges
%   2. contain no repeate
%   3. contain > 1 vertices
for i=1:length(F)
    % > 1 vxc
    if length(F{i}) < 2
        F{i} = [];
    elseif length(unique(F{i})) ~= length(F{i})
        F{i} = [];
    else
        Fe = char(join(string(F{i}), "-"));
        if sum(contains(T, Fe)) > 0
            F{i} = [];
        end
        %c = strfind(T,Fe);
        %if sum(cellfun(@isempty, c)) < m
        %    F{i} = [];
        %end
    end
end

% Remove empty false edges
F = F(find(cellfun(@isempty, F) == 0));

end

%% Test
%{
C = {char(join(string([1,2,3]), "-")), char(join(string([1,2]), "-")), char(join(string([3,4]), "-"))}

contains(C, string([3,4]))

find(C == [3,4])


cellfun(@contains, C, char(join(string([3,4]),"-"))

c = strfind(C,char(join(string([3,4]),"-")))
if cellfun(@isempty, c) < m
end

CC = c{:,:}

%}

