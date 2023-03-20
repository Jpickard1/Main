

function mouseHGparfor(tiArr)
% tiArr = 3

addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/'))
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/'))
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/'))

types = ["feed1", "fast", "feed2"];
type = types(tiArr);

D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');
[M, idxs] = HAT.multicorrelations(D, 3, 'Wang');
n = size(D,2);
% figure; histogram(M)

t = 0.8325; % Disconnected near 8355
hyperedges = idxs(find(M>t),:);
IM = HAT.hyperedge2IM(hyperedges);
HG = Hypergraph('IM', IM);
figure; plot(graph(HG.cliqueGraph)) % To check the graph is fully connected

% max(M)
% sum(M>0.85)


%%
rngs = [1 1508;
        1509 3017;
        3018 size(D,1)];

Ot = cell(size(D,1),1);
Dt = cell(size(D,1),1);
parfor t=rngs(tiArr, 1):rngs(tiArr,2)
    O = HGObsvNum(HG, D(t,:)');
    [d, ~] = greedyMON(O, n);
    Ot{t} = O;
    Dt{t} = d;
    % fileName = "mouseHG/" + string(type) + ".mat";
    % cmd = "save " + fileName + " " + "Ot Dt -v7.3";
    % eval(cmd);
    % disp(cmd);
end
fileName = "mouseHG/" + string(type) + "_par.mat";
cmd = "save " + fileName + " " + "Ot Dt -v7.3";
eval(cmd);
disp(cmd);



end
