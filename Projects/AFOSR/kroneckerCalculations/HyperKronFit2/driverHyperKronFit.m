function driverHyperKronFit(NameValueArgs)
%DRIVERKRONFIT Executes the KronFit algorithm on Great Lakes
%
%   PURPOSE: Execute the KronFit algorithm on graphs and hypergraphs
%       - Graphs: directed
%       - Hypergraphs: undirected
%
%   Arguments:
%       - system: GreatLakes (GL), lab computer (DBTM), laptop (LT)
%       - filePath: location of the file to read
%       - theta0: initiator matrix
%       - n0: size of initiator matrix
%       - firstPerm: number of initial permutations to try
%       - gradSamples: number of samples per gradient update
%       - learningRate: step size
%       - maxItrs: number of times to sample the gradient
%       - outputPath: file name/path of the results
%
%   driverKronFit('system','DBTM', 'n0', 2, 'filePath', 'C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt', 'firstPermItrs', 10000, 'gradSamples', '10000', 'maxItrs', '5');
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 19, 2023

%% Parse input

arguments
    NameValueArgs.system;
    NameValueArgs.filePath;
    NameValueArgs.theta0;
    NameValueArgs.n0;
    NameValueArgs.firstPermItrs;
    NameValueArgs.gradSamples;
    NameValueArgs.learningRate;
    NameValueArgs.maxItrs;
    NameValueArgs.outputPath;
end

% Required Inputs
sys           = NameValueArgs.system;
filePath      = NameValueArgs.filePath;

if isfield(NameValueArgs, 'theta0')
    theta0        = NameValueArgs.theta0;
    n0 = size(theta0,1);
elseif isfield(NameValueArgs, 'n0')
    n0            = NameValueArgs.n0;
    theta0 = rand(n0, n0);
else
    warning(['No defualt values or size is given for initial theta, ' ...
        'so a random initiator will be set with n0 = 2.']);
    n0 = 2;
    theta0 = rand(n0, n0);
end

if isfield(NameValueArgs, 'firstPermItrs')
    firstPermItrs = NameValueArgs.firstPermItrs;
else
    firstPermItrs = 1000;
end

if isfield(NameValueArgs, 'gradSamples')
    gradSamples = NameValueArgs.gradSamples;
else
    gradSamples = 50000;
end

if isfield(NameValueArgs, 'learningRate')
    learningRate = NameValueArgs.learningRate;
else
    learningRate = 1e-7;
end

if isfield(NameValueArgs, 'maxItrs')
    maxItrs = NameValueArgs.maxItrs;
else
    maxItrs = 10;
end

if isfield(NameValueArgs, 'outputPath')
    outputPath = NameValueArgs.outputPath;
else
    outputPath = filePath + "_HyperKronFit.out.mat";
end

% Variables still hardcodes
eps = 1e-4;

%% Preamble

if strcmp(sys, 'GL')
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/Code/utils'));
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit2/'));
elseif strcmp(sys, 'DBTM')
    addpath(genpath('C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations/HyperKronFit2/'));
    addpath(genpath('C:\Joshua\MissingData\Code\utils'));
elseif strcmp(sys, 'LT')
    addpath(genpath('C:\Users\picka\Documents\my_rojects\DBTM\Main\'));
else
    error(['Invalid System: the addpaths are configured for GreatLakes, the lab' ...
        'computer, or Joshuas laptop.']);
end

%% Execute Code
E = readAdjList(filePath, 0);       % Read file of adjacency list

if size(theta0, 1) ~= size(E,2)
    theta0 = rand(n0 * ones(size(E,2),1));
end

% Execute KronFit code
[theta, likelihoods, thetas] = HyperKronFit('E', E, 'theta0', theta0, 'maxItrs', maxItrs, 'gradSamples', gradSamples, 'firstPermItrs', firstPermItrs, 'learningRate', learningRate, 'eps', eps, 'v', true, 'p', false, 'debug', false);

%% Output
save(outputPath, "theta", "thetas", "likelihoods");

end


% NaiveKronFitArgs = struct;
% NaiveKronFitArgs.A = A;
% NaiveKronFitArgs.theta0 = theta0;
% NaiveKronFitArgs.maxItrs = maxItrs;
% NaiveKronFitArgs.gradSamples = gradSamples;
% NaiveKronFitArgs.firstPermItrs = firstPermItrs;
% NaiveKronFitArgs.learningRate = learningRate;
% NaiveKronFitArgs.eps = eps;
% NaiveKronFitArgs.v = true;
% NaiveKronFitArgs.p = false;
% NaiveKronFitArgs.debug = false;