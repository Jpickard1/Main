% Equation 13
function similarity=katz_index(adj, beta)
    largest_eigenvalue = eigs(adj, 1);
    if beta > (1/largest_eigenvalue)
        disp('ERROR: Beta is too large. Beta must be < 1/largest eigenvalue of adj.');
        disp(str(beta) + '>' + str(1/largest_eigenvalue));
    end
    similarity = inv(eye(size(adj)) - (beta * adj)) - eye(size(adj));
end
