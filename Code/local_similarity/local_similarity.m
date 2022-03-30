% Description: This function computes a similarity matrix between each
% vertex in a network. It accepts and adjacency matrix as well as a
% parameter for which local similarity score should be used. These scores
% are defined in the paper: "Uncovering Missing Links with Cold Ends" and
% currently available ones include:
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
% See: https://drive.google.com/file/d/1TsXGWMqqt5euCtxnEQBodo3FM74gyvzO/view
function similarity=local_similarity(adj, index)
    similarity = zeros(size(adj));
    for i=1:height(adj)
        for j=i+1:height(adj)
            % i_j_similarity = 0;
            switch index
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
            end
            similarity(i, j) = i_j_similarity;
            similarity(j, i) = i_j_similarity;
        end
    end
end
