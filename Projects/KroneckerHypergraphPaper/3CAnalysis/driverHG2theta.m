function driverHG2theta(NameValueArgs)
%DRIVERKRONFIT Executes the KronFit algorithm on Great Lakes
%
%   PURPOSE: Execute the KronFit algorithm on graphs and hypergraphs
%
%   Arguments:
%       - worker: array job number
%       - system: GreatLakes (GL), lab computer (DBTM), laptop (LT)
%       - filePath: location of the file to read
%       - n0: size of initiator tensor
%       - firstPerm: number of initial permutations to try
%       - gradSamples: number of samples per gradient update
%       - learningRate: step size
%       - maxItrs: number of times to sample the gradient
%       - outputPath: file name/path of the results
%       - epsilon: epsilon value used in hypergraph construction (for file
%       selection)
%
%   driverHG2theta('system','DBTM', 'n0', 2, 'filePath', 'C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt', 'firstPermItrs', 10000, 'gradSamples', '10000', 'maxItrs', '5');
%
%   filePath = 'C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt'
%
%   filePath = "/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/";
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 21, 2023

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

disp('EXECUTE: driverHipHop2HG(sys, worker, epsilon, path2data, path2out)')
disp(['         WARNING: hyperedges must contain vertices at least 7 spots' ...
      ' appart on the polymer to be considered. This is hardcoded']);
disp('         input arguments are printed out in order:');
disp(NameValueArgs);

% Required Inputs
sys           = NameValueArgs.system;
path2data      = NameValueArgs.filePath;

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
    path2out = NameValueArgs.outputPath;
else
    error('Specify the output location');
end

% Variables still hardcoded
eps = 1e-4;

%% Preamble

if strcmp(sys, 'GL')
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit2/'));
elseif strcmp(sys, 'DBTM')
    addpath(genpath('C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations/HyperKronFit2/'));
elseif strcmp(sys, 'LT')
    addpath(genpath('C:\Users\picka\Documents\my_rojects\DBTM\Main\Projects\AFOSR\kroneckerCalculations\HyperKronFit2\'));
else
    error(['Invalid System: the addpaths are configured for GreatLakes, the lab' ...
        'computer, or Joshuas laptop.']);
end

%% Execute Script

% 1. Get File Path
path2data = path2data + "adjList_";
if worker <= 4
    mid = "OFF";
elseif worker <= 9
    mid = "ON";
else
    mid = "HIGH";
end
path2data = path2data + mid + "_";
path2out  = path2out  + "HKF_" + mid + "_";
suffixTXT = ".txt";
suffixMAT = ".mat";

% 2. Set file range
fileRange = mod(worker, 5) * 40 + [1:40];

% 3. Iterate over all files
for f=1:length(fileRange)
    % Set file path
    filePath = path2data + string(fileRange(f)) + "_" + string(epsilon) + suffixTXT;
    outPath = path2out + string(fileRange(f)) + "_" + string(epsilon) + suffixMAT;
    disp(filePath);

    % Read file of adjacency list
    E = readAdjList(filePath, 0);
    
    % Filter out hyperedges that do not contain meaningful contact information
    % regarding the structure of the polymer. While subjective, we define an
    % interesting hyperedge as containing beads that are further than 7 beads
    % away from one another.
    minDistanceOfInterest = 7;
    keep = find(range(E,2) >= minDistanceOfInterest);
    E = E(keep, :);
    
    % Execute KronFit code
    [theta, likelihoods, thetas] = HyperKronFit('E', E, 'theta0', theta0, 'maxItrs', maxItrs, 'gradSamples', gradSamples, 'firstPermItrs', firstPermItrs, 'learningRate', learningRate, 'eps', eps, 'v', true, 'p', false, 'debug', false);
    
    % Output
    save(outPath, "theta", "thetas", "likelihoods");

end


end