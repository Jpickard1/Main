% Description: This function takes in edges and partitions them according
% to various schemes. Current partitioning schemes are:
%   - 'rand': randomly splits edges between training and testing
%   - 'pop': splits them into multiple groups according to popularity
function [known, unknown]=popularity_partition_edges(popularity_matrix, D)
    adj = logical(popularity_matrix);
    popularities = sort(popularity_matrix(:));
    noedge = find(popularities == 0);
    popularities = popularities(noedge(end) + 1:end);
    edges_per_group = round(length(popularities) / D);
    % Make cell arrays containing known and unknown edges for each
    % popularity group
    known = cell(D,1);
    unknown = cell(D,2);
    % Iterate over popularity levels
    for group=1:D
        pop_low = popularities(max([edges_per_group * (group - 1), 1]));
        pop_high = popularities(min([edges_per_group * group, length(popularities)]));
        %if group == D
        %    pop_high = popularities(end);
        %end
        % Set the known and unknownedges
        pop_unknown = (popularity_matrix < pop_high) & (popularity_matrix > pop_low);
        pop_known = (adj & ~pop_unknown);
        % Save them in the network
        known{group} = pop_known;
        unknown{group} = pop_unknown;
    end
end
