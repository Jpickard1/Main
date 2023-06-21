function scratchEnron(sys)


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


E = readAdjList("adjLists/email-Enron_3.txt");
n = max(max(E));

D = zeros(n,1);
for i=1:n
    D(i) = sum(E == i, 'all');
end

figure; plot(sort(D));

theta0 = rand(3,3,3);
[theta, likelihoods, thetas] = HyperKronFit('E',E,'firstPermItrs',10000,'gradSamples',100000,'maxItrs',100,'theta0',theta0,'learningRate',1e-7);

path2out = "enron-scratch-3-3.mat";
save(path2out, "theta", "likelihoods", "thetas")

end