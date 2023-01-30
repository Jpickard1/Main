function similarity=compute_similarity(adj, index)
% Description: This function computes a similarity matrix between each
% vertex in a network. It accepts and adjacency matrix as well as a
% parameter for which similarity index should be used.
%
% Local Similarity Measures:
%  - AA: Adamic Adar
%  - CN: Common Neighbors
%  - HD: Hub Depressed
%  - HP: Hub Promoted
%  - JC: Jaccard
%  - LHN: Leicht-Holme-Newman
%  - PA: Preferential Attachment
%  - RS: Resource Allocation
%  - SA: Salton
%  - SO: Sorensen
% Global Similarity Measures:
%  - 'K': Katz
%  - 'ACT': Average Commute Time
%  - 'RWR': Random Walk with Restart
%
% References:
% Local indices: "Uncovering Missing Links with Cold Ends"
%    https://drive.google.com/file/d/1TsXGWMqqt5euCtxnEQBodo3FM74gyvzO/view
% Global indices: "Link prediction in complex networks: A survey"
%    https://drive.google.com/file/d/1ORRPl2rUiboOJ5FdYnaIvl3uk4lAT2Ff/view
    similarity = zeros(size(adj));
    if strcmp(similarity, 'K')
        similarity = katz_index(adj);
    end
    jIDXS = 1:height(adj);
    parfor i=1:height(adj)
        i_j_similarity = -1;
        for j=jIDXS
            if j < i+1
                continue
            end
            % i_j_similarity = 0;
            switch index
                % Local similarity
                case 'AA'
                    i_j_similarity = adamicAdar_index(adj, i, j);
                case 'CN'
                    i_j_similarity = commonNeighbors_index(adj, i, j);
                case 'HD'
                    i_j_similarity = hubDepressed_index(adj, i, j);
                case 'HP'
                    i_j_similarity = hubPromoted_index(adj, i, j);
                case 'JC'
                    i_j_similarity = jaccard_index(adj, i, j);
                case 'LHN'
                    i_j_similarity = leichtHolmeNewman_index(adj, i, j);
                case 'PA'
                    i_j_similarity = preferentialAttachment(adj, i, j);
                case 'RS'
                    i_j_similarity = resourceAllocation(adj, i, j);
                case 'SA'
                    i_j_similarity = salton_index(adj, i, j);
                case 'SO'
                    i_j_similarity = sorensen_index(adj, i,j);
                % Global similarity
                case 'ACT'
                    i_j_similarity = averageCommuteTime_index(adj, i, j);
                case 'RWR'
                    i_j_similarity = randomWalkWithRestart_index(adj, i, j);
                % case 'K'
                    % similarity = katz_index(adj);
                    % i = 2 * height(adj);
                    % j = i;
                    % return
            end
            % similarity(i, j) = i_j_similarity;
            similarity(j, i) = i_j_similarity;
        end
    end
    similarity = similarity + similarity';
end
