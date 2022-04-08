function scree(matrix)
% SCREE: This function generates a scree plot for square matrices
%
%   Joshua Pickard jpic@umich.edu
%   April 8, 2022
    eigenvalues = real(eig(matrix));
    plot(eigenvalues, '-o');
    ylabel('Eigenvalue');
    xlabel('Component');
    title('Scree Plot');
end
