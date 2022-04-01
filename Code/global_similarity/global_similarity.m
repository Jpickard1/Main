% Description: This function computes a similarity matrix between each
% vertex in a network. It accepts and adjacency matrix as well as a
% parameter for which local similarity score should be used. These scores
% are defined in the paper: "Link prediction in complex networks: A survey"
% and currently available ones include:
%  - 'K': Katz
%  - 'ACT': Average Commute Time
%  - 'RWR': Random Walk with Restart
% See: https://drive.google.com/file/d/1ORRPl2rUiboOJ5FdYnaIvl3uk4lAT2Ff/view
function similarity=global_similarity(adj, index)
    if strcmp(index, 'K')
        similarity = katz_index(adj);
    else
        similarity = zeros(size(adj));
        for i=1:height(adj)
            for j=i+1:height(adj)
                % i_j_similarity = 0;
                switch index
                    case 'ACT'
                        i_j_similarity = averageCommuteTime_index(adj, i, j);
                    case 'RWR'
                        i_j_similarity = randomWalkWithRestart_index(adj, i, j);
                end
                similarity(i, j) = i_j_similarity;
                similarity(j, i) = i_j_similarity;
            end
        end
    end
end
