function matrix= discretize_scHiCluster(matrix)
    %% Set parameters
    p_top = 0.8; % The top percent of links that should be included
    
    %% Select threshhold value
    matrix_values = sort(matrix(:));
    threshold = matrix_values(round(p_top * length(matrix_values)));
    
    %% Discretize
    matrix = (matrix > threshold);
end
