addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/'))
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/'))
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/'))

D = csvread('Copy_of_mouse 44_fed_fasted_refed traces.csv');
n = size(D, 2);

rngs = [1 1508;
        1509 3017;
        3018 size(D,1)];

% neurons = ones(15, 1); % neurons([8, 10, 13],:) = 0;
neurons = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 0 1 0 0 1 0];
neurons = logical(neurons);
D1 = D(rngs(1,1):rngs(1,2), neurons);
D2 = D(rngs(2,1):rngs(2,2), neurons);
D3 = D(rngs(3,1):rngs(3,2), neurons);
[M1, idxs] = HAT.multicorrelations(D1, 3, 'Wang');
[M2, ~] = HAT.multicorrelations(D2, 3, 'Wang');
[M3, ~] = HAT.multicorrelations(D3, 3, 'Wang');
figure; hold on; histogram(M1, 20); histogram(M2, 20); histogram(M3, 20);
t = 0.95; % Disconnected near 8355
hyperedges1 = idxs(find(M1>t),:);
IM1 = HAT.hyperedge2IM(hyperedges1);
HG1 = Hypergraph('IM', IM1);
hyperedges2 = idxs(find(M2>t),:);
IM2 = HAT.hyperedge2IM(hyperedges2);
HG2 = Hypergraph('IM', IM2);
hyperedges3 = idxs(find(M3>t),:);
IM3 = HAT.hyperedge2IM(hyperedges3);
HG3 = Hypergraph('IM', IM3);
figure; plot(graph(full(HG1.cliqueGraph)));
figure; plot(graph(full(HG3.cliqueGraph)));
figure; plot(graph(full(HG2.cliqueGraph)));

HG = {HG1, HG2, HG3}
save mouseHGs2.mat HG

%% Linear MON identification

for hgNum=1:3
    load('mouseHGs2.mat');
    HG = HG{hgNum};

    % Reduce HG to only connected components
    degree = sum(full(HG.IM), 2);
    vxc = find(degree > 0);
    IM = full(HG.IM);
    IM = IM(vxc,:);
    HG = Hypergraph('IM', IM);

    A = full(HG.cliqueGraph);

    O = cell(size(A,1), 1);
    for vx=1:size(A,1)
        C = zeros(1, size(A,1)); C(vx) = 1;
        O{vx} = obsv(A, C);
    end    
    [D, ~] = greedyMON(O, length(vxc));             % Greedy Node Selection
    disp(length(O))
    disp(D);
    
end

%%

% Compute observability matrices for all vertices
O = cell(n,1);                  % Cell to hold observability matrices for all vertices
symVars = symvar(sym('x_%d',[n 1]));  % Get symbolic variables 
for vx=1:n
    disp(vx);
    Ci = zeros(1,n); Ci(vx) = 1;% Equation 9
    Oi = cell(7,1);             % Compute first equality in equation 10
    for i=1:7
        Oi{i} = Ci * J{i};      % Compute first equality in equation 10
    end
    Oimat = sym([7,n]);         
    for i=1:7                   % Set symbolic matrix to save equation 10 for vertex vx
        for j=1:n               % Compute second equality in equation 10
            Oimat(i,j) = gradient(Oi{i}, symVars(j));
        end
    end
    O{vx} = Oimat;              % Save observability matrix for specific vertex
end

%% Aligning the vertices in matlab with those for the latex plots
% resetdefaultpath
clear
load('mouseHGs2.mat')
h = HG{2};
sum(full(h.cliqueGraph))

h = HG{1};
sum(full(h.cliqueGraph))