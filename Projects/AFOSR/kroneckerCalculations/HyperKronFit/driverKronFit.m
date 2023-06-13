function driverKronFit(NameValueArgs)
%DRIVERKRONFIT Executes the KronFit algorithm on Great Lakes
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
%
%   driverKronFit('system','DBTM', 'n0', 2, 'filePath', 'C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations\HyperKronFit\as20graph.txt', 'firstPermItrs', 10000, 'gradSamples', '10000', 'maxItrs', '5');
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 13, 2023

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
        'so a random initiator is set with n0 = 2.']);
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

% Variables still hardcodes
eps = 1e-3;

%% Preamble

if strcmp(sys, 'GL')
    addpath(genpath(''));
elseif strcmp(sys, 'DBTM')
    addpath(genpath('C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations'));
    addpath(genpath('C:\Joshua\MissingData\Code\utils'));
elseif strcmp(sys, 'LT')
    addpath(genpath('C:\Users\picka\Documents\my_rojects\DBTM\Main\'));
else
    error(['Invalid System: the addpaths are configured for GreatLakes, the lab' ...
        'computer, or Joshuas laptop.']);
end

%% Execute Code
E = readAdjList(filePath, 0);       % Read file of adjacency list
A = sparse(E(:,1), E(:,2), 1);      % Convert adjacency list to adjacency
                                    % matrix

% Execute KronFit code
[theta, likelihoods, thetas] = NaiveKronFit('A', A, 'theta0', theta0, ...
    'maxItrs', maxItrs, 'gradSamples', gradSamples, 'firstPermItrs', ...
    firstPermItrs, 'learningRate', learningRate, 'eps', eps, 'v', true, ...
    'p', false, 'debug', false);
% [theta, likelihoods] = NaiveKronFit(A, true, true, 2, theta0, 10);  

%% Output


end