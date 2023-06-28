function A_avg = KronFilter4(A, k, blockSizes, numIterations)
% Perform high-rank Kronecker decomposition using binning and block averaging
% Inputs:
%   A: Original matrix
%   k: Desired rank
%   blockSizes: Vector specifying the sizes of the blocks in each dimension
%   numIterations: Number of iterations for refinement
% Output:
%   A_avg: High-rank Kronecker decomposition approximation

% Get the dimensions of the original matrix
[m, n] = size(A);

% Initialize A_avg as the original matrix
A_avg = A;

% Perform block averaging
for i = 1:numel(blockSizes)
    blockSize = blockSizes(i);
    numBlocks = floor(m/blockSize);
    
    % Reshape matrix into blocks and compute block averages
    A_avg = blockproc(A_avg, [blockSize blockSize], @(x) mean(x.data(:)), 'PadPartialBlocks', true);
end

% Repeat refinement steps until the desired rank is achieved
for iter = 1:numIterations
    % Generate a low-rank approximation matrix
    B = randn(k);
    
    % Compute Kronecker approximation
    A_new = kron(A_avg, B);
    
    % Refinement: Update B using least squares
    B = lsqnonneg(kron(eye(k), A_avg), A(:));
    
    % Reshape A_new to match the original matrix size
    A_new = reshape(A_new, size(A));
    
    % Update A_avg
    A_avg = A_new;
end

end
