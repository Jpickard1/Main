% Description: This function predicts missing edges between the most
% similar nodes that do not yet have an edge between them. The similarity
% matrix should be directly output by the funciton compute_similarity.
% Additionally, a parameter should be set to specify the number of missing
% links expected, so that only the top n% of similar links are predicted.

function predicted_network = similarity_predict(known_network, similarity, n_percent)
    % Create predicted network to contain all known links
    predicted_network = false(size(known_network));

    % Set similarity between known neighbors to 0 so that they are not
    % predicted
    similarity(known_network == true) = 0;

    % Deterime the number of known edges to predict
    n_known_edges = sum(sum(known_network));
    n_edges_to_predict = round(n_percent * n_known_edges);

    min_similarity = min(min(similarity));

    % Predict Edges
    e = 1;
    while e < n_edges_to_predict
        % Find maximum similarity value
        m = max(max(similarity));
        
        % Find the nodes this similarity is between
        [i,j] = find(similarity == m);
        edges = [i j];

        % Set predict the link
        for edge=1:height(edges)
            predicted_network(edges(edge,1), edges(edge,2)) = true;
        end
        % predicted_network(i,j) = true;

        % Set similarity to minimum so that it is not used again
        similarity(i,j) = min_similarity;

        e = e + length(i);
    end
    predicted_network = (predicted_network | known_network);
end
