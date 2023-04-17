function scratch2(hgNum)
    disp('Start');
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Main/'))
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/'))
    addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/'))
    load("mouseHGs2.mat")

    HG = HG{hgNum};

    % Reduce HG to remove disconnected vertices
    degree = sum(full(HG.IM), 2);
    vxc = find(degree > 0);
    IM = full(HG.IM);
    IM = IM(vxc,:);
    HG = Hypergraph('IM', IM);

    O = HGObsvSym(HG);
    disp('Observability Calculations Done');
    if hgNum == 1
        save O1.mat O;
    elseif hgNum == 2
        save O2.mat O;
    elseif hgNum == 3
        save O3.mat O;
    end
    [D, ~] = greedyMON(O, length(vxc));             % Greedy Node Selection
    disp(D);
    
end
