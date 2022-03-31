function matrix=linearConvolution(matrix, kernel_size)
    if nargin == 1
        kernel_size = 3;
    end
    kernel = ones(kernel_size);
    kernel = kernel / sum(sum(kernel));
    output = conv2(matrix, kernel, 'same')
end