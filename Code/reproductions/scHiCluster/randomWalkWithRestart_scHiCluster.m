function Q = randomWalkWithRestart_scHiCluster(matrix)
    % Parameters
    restart_probability = 0.05;
    tol = 1e-6;

    % Row normalize the matrix
    % matrix = matrix / sum(matrix,2);
    matrix = normr(matrix);
    
    % Random Walk With Restart
    Q = eye(size(matrix));
    matrix_norm = tol + 1;
    while matrix_norm > tol
        Q_new = (1-restart_probability)*Q*matrix + restart_probability * eye(size(matrix));
        matrix_norm = norm(Q_new - Q);
        Q = Q_new;
    end
end
