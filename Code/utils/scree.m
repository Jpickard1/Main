function scree(matrix)
% SCREE: This function generates a scree plot for square matrices
%
%   Joshua Pickard jpic@umich.edu
%   April 8, 2022
    % Row normalize the matrix
    for i=1:length(matrix)
        if sum(matrix(i,:)) ~= 0
            matrix(i,:) = matrix(i,:) / sum(matrix(i,:));
        end
    end
    eigenvalues = flip(sort(real(eig(matrix))));
    plot(eigenvalues, '-o');
    ylabel('Eigenvalue');
    xlabel('Component');
end
