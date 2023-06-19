function [theta, likelihoods, thetas] = HyperKronFit(NameValueArgs)
%HYPERKRONFIT Summary of this function goes here
%   Detailed explanation goes here
%{
    filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt";
    E =  readAdjList(filePath, 0);       % Read file of adjacency list
    NameValueArgs = struct;
    NameValueArgs.E             =  readAdjList(filePath, 0);       % Read file of adjacency list
    NameValueArgs.theta0        = rand(2,2);
    NameValueArgs.maxItrs       = 25;
    NameValueArgs.gradSamples   = 50000;
    NameValueArgs.firstPermItrs = 10000;
    NameValueArgs.learningRate  = 1e-6;
    NameValueArgs.eps           = 1e-3;
    NameValueArgs.debug         = false;
    NameValueArgs.v             = true;
    NameValueArgs.p             = false;
    HyperKronFit(NameValueArgs)
    HyperKronFit("E",E,"theta0",rand(2,2),"maxItrs",25,"firstPermItrs",10000,"learningRate",1e-5,"gradSamples",50000,"v",true)
%}

%% Parse parameters
arguments
    NameValueArgs.E;
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

E = NameValueArgs.E;
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
n  = max(max(E));
k  = size(E,2);

likelihoods = zeros(maxItrs, 1);
for itr=1:maxItrs
    if verbose
        fprintf("Itr: %d ]\n", itr);
    end

    % Evaluate likelihood and gradient
    [l, gradients] = sampleGradient(theta, gradSamples, firstPermItrs, E);
    % Update model parameters
    thetaOld = theta;
    theta = theta + learningRate * gradients;

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

