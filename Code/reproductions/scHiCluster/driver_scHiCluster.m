function matrix = driver_scHiCluster(matrix)
    matrix = linearConvolution_scHiCluster(matrix);
    matrix = randomWalkWithRestart_scHiCluster(matrix);
    matrix = discretize_matrix(matrix);
end