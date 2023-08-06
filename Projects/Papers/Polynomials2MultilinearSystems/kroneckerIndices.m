function indices = kroneckerIndices(n, d, i)
    indices = zeros(1, d);

    for j = d:-1:1
        indices(d - j + 1) = mod(floor((i - 1) / n^(j - 1)), n) + 1;
    end
end
