function [theta, likelihoods, thetas] = NaiveKronFit(NameValueArgs)
%NAIVEKRONFIT This function is slow but does a brute force evaluation of
% the kron fit problem according to equation 5.5 in Jure Leskovec's thesis.
%
%   Arguments:
%       - A: adjacency matrix
%       - theta0: initial guess for kronecker parameters
%       - maxItrs: maximum number of iterations
%       - gradSamples: number of samples per gradient approximation
%       - firstPermItrs: number of permutations to try initially
%       - learningRate (lr): learning rate for each gradient step
%       - eps: minimum distance of values in theta from 0 and 1 (default
%       1e-2)
%       - verbose (v): flag to surpress or generate output (default is true)
%       - debug: flag for debugging
%       - plotting (p): flag to plot likelihood function every 50
%       iterations
%
%   TODO: Improvements/bug fixes that need to be made
%       - modify the code to work for hypergraphs (Joshua)
%       - find error in gradient calculation which leads to the drift from
%       the correct parameters to 1 (Joshua)
%       - change how the learning rate is set (Vivan)
%       - design a plan for when this algorithm/code should terminate
%       (Vivan)
%       - modify the code to periodically save the results (Vivan)
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 5, 2023

%% Parse parameters
arguments
    NameValueArgs.A;
    NameValueArgs.theta0;
    NameValueArgs.maxItrs;
    NameValueArgs.gradSamples;
    NameValueArgs.firstPermItrs;
    NameValueArgs.learningRate;
    NameValueArgs.eps;
    NameValueArgs.debug;
    NameValueArgs.v;
    NameValueArgs.verbose;
    NameValueArgs.p;
end
A = NameValueArgs.A;
theta = NameValueArgs.theta0;
if isfield(NameValueArgs, 'maxItrs')
    maxItrs = NameValueArgs.maxItrs;
else
    maxItrs = 5;
end
if isfield(NameValueArgs, 'gradSamples')
    gradSamples = NameValueArgs.gradSamples;
else
    gradSamples = 100000;
end
if isfield(NameValueArgs, 'firstPermItrs')
    firstPermItrs = NameValueArgs.firstPermItrs;
else
    firstPermItrs = 10000;
end
if isfield(NameValueArgs, 'learningRate')
    learningRate = NameValueArgs.learningRate;
elseif isfield(NameValueArgs, 'lr')
    learningRate = NameValueArgs.lr;
else
    learningRate = 1e-7;
end
if isfield(NameValueArgs, 'eps')
    eps = NameValueArgs.eps;
else
    eps = 1e-2;
end
if isfield(NameValueArgs, 'debug')
    debug = NameValueArgs.debug;
else
    debug = false;
end
if isfield(NameValueArgs, 'verbose')
    verbose = NameValueArgs.verbose;
elseif isfield(NameValueArgs, 'v')
    verbose = NameValueArgs.v;
else
    verbose = false;
end
if isfield(NameValueArgs, 'p')
    plotOn = NameValueArgs.p;
else
    plotOn = false;
end


% Output arguments
if nargout == 3
    thetas = cell(maxItrs + 1, 1);
    thetas{1} = theta;
end

%% Execute code
n0 = size(theta, 1);
n  = size(A,1);
k  = length(size(A));
E = getEdgesFromAdj(A);

directed = false; % This will remain hardcoded for now
if k == 2; directed = true; end

likelihoods = zeros(maxItrs, 1);
for itr=1:maxItrs
    if verbose
        fprintf("Itr: %d ]\n", itr);
    end
    
    % Evaluate likelihood and gradient
    % [l, gradients] = evaluateGradient(A, theta, debug);
    [l, gradients] = sampleGradient(A, theta, gradSamples, firstPermItrs, debug, directed, E);
    % Update model parameters
    thetaOld = theta;
    theta = theta + learningRate * gradients;

    % lr = 0.95 * lr;
    for i=1:n0
        for j=1:n0
            if theta(i,j) > 1 - eps; theta(i,j) = 1 - eps; end
            if theta(i,j) < eps
                theta(i,j) = eps;
            end
        end
    end

    if verbose
        fprintf("CurrentLL: %d\n", l);
        fprintf("Gradient Updates:\n");
        for i=1:numel(theta)
            fprintf("    %d]  %f = %f + %f  \t Grad:  %f\n", i, theta(i), thetaOld(i), learningRate * gradients(i), gradients(i));
        end
        fprintf(' \n');
    end

    % Outputs
    if plotOn && mod(itr, 50) == 0
        disp(itr);
        figure; plot(real(likelihoods));
        pause(1);
    end
    likelihoods(itr) = l;
    if nargout == 3
        thetas{itr + 1} = theta;
    end
end

end

