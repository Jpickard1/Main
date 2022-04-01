function matrix = scHiCluster_imputation(matrix)
    matrix = linearConvolution_scHiCluster(matrix);
    matrix = randomWalkWithRestart_scHiCluster(matrix);
    
end
