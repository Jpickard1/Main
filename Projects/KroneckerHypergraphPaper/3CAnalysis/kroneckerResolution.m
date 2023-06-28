%% Kronecker Filter to Hi-C data
clear
j = 1;
s = 1;
SIMPARMS = ["OFF","ON","HIGH"];
path2data = 'C:\Joshua\MissingData\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

figure; imagesc(D)

%%

[IMp, F] = kronFilter(D, [50 50])

figure;
imagesc(IMp)

figure;
imagesc(F{1})

figure;
imagesc(kron(kron(IMp, F{1}), F{2}))

%%
B = zeros(3,3);
figure;
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
path2data = 'C:\Joshua\MissingData\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
ctr = 1;
for s=1:3
    path2sims = path2data + SIMPARMS(s) + "\conf.";
    for j=1:nsamps
        filePath = path2sims + string(j) + ".DNA";
        T = readLAMMPSoutput(filePath);
        cords = getCords(T);
        subplot(3,nsamps,ctr); ctr = ctr + 1;
        for i = 1:size(B,1)
            vxc = (i-1) * 1000 + [1:1000];
            plot3(cords(vxc,1), cords(vxc,2), cords(vxc,3), '.-'); hold on;
        end
        title(SIMPARMS(s) + ": " + string(j));
    end
end

%%
D = squareform(pdist(cords));


figure;
subplot(1,3,1); imagesc(D)
subplot(1,3,2); imagesc(F1)
subplot(1,3,3); imagesc(D2)

nBins = 5;
sBins = size(D,1) / nBins;

F1 = zeros(nBins, nBins);
for i=1:nBins
    for j=1:nBins
        F1(i,j) = sum(D((i-1)*sBins+1:i*sBins,(j-1)*sBins+1:j*sBins), 'all');
    end
end

D2 = D;
for i=1:nBins
    for j=1:nBins
        D2((i-1)*sBins+1:i*sBins,(j-1)*sBins+1:j*sBins) = D((i-1)*sBins+1:i*sBins,(j-1)*sBins+1:j*sBins) ./ F1(i,j);
    end
end

F2 = zeros(nBins, nBins);
for i=1:nBins
    for j=1:nBins
        F2(i,j) = sum(D2((i-1)*sBins+1:i*sBins,(j-1)*sBins+1:j*sBins), 'all');
    end
end


%%
clearvars -except D

nBins = 5;
sBins = size(D,1) / nBins;

F1 = zeros(nBins, nBins);
B1 = cell(nBins, nBins);
for i=1:nBins
    for j=1:nBins
        F1(i,j) = sum(D((i-1)*sBins+1:i*sBins,(j-1)*sBins+1:j*sBins), 'all');
        B1{i,j} = D((i-1)*sBins+1:i*sBins,(j-1)*sBins+1:j*sBins);
    end
end

dist1 = zeros(nBins^2, nBins^2);
for i=1:nBins^2
    [x1, y1] = ind2sub([nBins, nBins], i);
    for j=1:nBins^2
        [x2, y2] = ind2sub([nBins, nBins], j);
        dist1(i,j) = norm(B1{x1,y1} - B1{x2,y2});
    end
end
figure; imagesc(dist1);

F2 = cell(nBins, nBins);
for i1=1:nBins
    for j1=1:nBins
        A = B1{i1,j1};
        sBins = size(A, 1) / nBins;
        F2ij = zeros(nBins, nBins);
        for i2=1:nBins
            for j2=1:nBins
                F2ij(i2,j2) = sum(A((i2-1)*sBins+1:i2*sBins,(j2-1)*sBins+1:j2*sBins), 'all');
            end
        end
        F2{i1,j1} = F2ij / sum(F2ij, 'all');
    end
end

figure; tiledlayout(5,5)
for i=1:nBins
    for j=1:nBins
        nexttile;
        imagesc(F2{i,j}); title(string(i) + " - " + string(j));
    end
end


%%
clearvars -except D

D1 = zeros(size(D,1) * 0.5, size(D,1) * 0.5);
P1 = zeros(2,2);
for i=1:size(D1,1)
    for j=1:size(D1,1)
        D1(i,j) = mean(D((i-1)*2+1:i*2,(j-1)*2+1:j*2), 'all');
        P1 = P1 + D((i-1)*2+1:i*2,(j-1)*2+1:j*2);
    end
end
P1 = P1 / (i*j);
D2 = zeros(size(D1,1) * 0.5, size(D1,1) * 0.5);
P2 = zeros(2,2);
for i=1:size(D2,1)
    for j=1:size(D2,1)
        D2(i,j) = mean(D1((i-1)*2+1:i*2,(j-1)*2+1:j*2), 'all');
        P2 = P2 + D1((i-1)*2+1:i*2,(j-1)*2+1:j*2);
    end
end
P2 = P2 / (i*j);
D3 = zeros(size(D2,1) * 0.5, size(D2,1) * 0.5);
P3 = zeros(2,2);
for i=1:size(D3,1)
    for j=1:size(D3,1)
        D3(i,j) = mean(D2((i-1)*2+1:i*2,(j-1)*2+1:j*2), 'all');
        P3 = P3 + D2((i-1)*2+1:i*2,(j-1)*2+1:j*2);
    end
end
P3 = P3 / (i*j);
D4 = zeros(floor(size(D3,1) * 0.5), floor(size(D3,1) * 0.5));
P4 = zeros(2,2);
for i=1:size(D4,1)
    for j=1:size(D4,1)
        D4(i,j) = mean(D3((i-1)*2+1:i*2,(j-1)*2+1:j*2), 'all');
        P4 = P4 + D3((i-1)*2+1:i*2,(j-1)*2+1:j*2);
    end
end
P4 = P4 / (i*j);

% Reconstruct image
R4 = kron(D4, P4);
R3 = kron(R4, P3);
R2 = kron(R3, P2);
R1 = kron(R2, P1);

figure; layers = 5;
subplot(2,layers,1); imagesc(D);  title("Level: 0");
subplot(2,layers,2); imagesc(D1); title("Level: 1");
subplot(2,layers,3); imagesc(D2); title("Level: 2");
subplot(2,layers,4); imagesc(D3); title("Level: 3");
subplot(2,layers,5); imagesc(D4); title("Level: 4");
subplot(2,layers,6);  imagesc(R4); title("Level: 4");
subplot(2,layers,7);  imagesc(R3); title("Level: 3");
subplot(2,layers,8);  imagesc(R2); title("Level: 2");
subplot(2,layers,9);  imagesc(R1); title("Level: 1");

figure; layers = 4;
subplot(1,layers,1); imagesc(P1);
subplot(1,layers,2); imagesc(P2);
subplot(1,layers,3); imagesc(P3);
subplot(1,layers,4); imagesc(P4);



%% Make some example images with Kronecker
clear
A = sin(0:0.5:2*pi);
A = repmat(A,[length(A), 1]);
figure; imagesc(A)

B = kron(A, A);
figure; imagesc(B);

%%
[IMp, F] = kronFilter(B,13);

figure;
subplot(1,2,1); imagesc(IMp)
subplot(1,2,2); imagesc(A)

Br = kron(IMp, F{1});
Br = (max(max(B)) / max(max(Br))) * Br;
delta = B - Br;

figure;
subplot(1,3,1); imagesc(B);
subplot(1,3,2); imagesc(Br);
subplot(1,3,3); imagesc(delta);


norm(delta)
norm(B)


%%

[IMp, F] = kronFilter(D, 2*ones(10,1));
% [IMp, F] = kronFilter(rand(2^5, 2^5), 2*ones(4,1));

function [IMp, F] = kronFilter(IM, levelSizes)

    levelSize = levelSizes(1);

    IMp = zeros(floor(size(IM,1) * (1/levelSize)), floor(size(IM,1) * (1/levelSize)));
    F = zeros(levelSize, levelSize);
    for i=1:size(IMp,1)
        for j=1:size(IMp,1)
            IMp(i,j) = mean(IM((i-1)*levelSize+1:i*levelSize,(j-1)*levelSize+1:j*levelSize), 'all');
            F = F + IM((i-1)*levelSize+1:i*levelSize,(j-1)*levelSize+1:j*levelSize);
        end
    end
    F = {F / (i*j)};

    if length(levelSizes) > 1
        [IMp, filters] = kronFilter(IMp, levelSizes(2:end));
        FF = cell(length(filters) + 1,1);
        FF{length(FF)} = F{1};
        % if length(filters) > 1
            for i=1:length(filters)
                FF{i} = filters{i};
            end
        % else
        %     FF{1} = filters;
        % end
        F = FF;
    end

end

function IM = kronExpansion(IMp, filters)
    
    IM = IMp;
    for i=1:lenght(filters)
        IM = kron(IM, filters{i});
    end

end



















