% load a matrix D

clearvars -except D
[U, S, V] = svd(D);
sigmas = diag(S);
sigmas(1:10)

[IM, F] = KronFilter(D, 2*ones(10,1));


%%
A = rand(2^5,2^5);
[IM, F] = KronFilter(A, [2 2 2]);

%%
figure; hold on;
for i=1:10
    [u,s,v] = svd(IM{i});
    sigmas = diag(s);
    plot(sigmas(1:min([length(sigmas), 20])));
end
[U, S, V] = svd(D);
sigmas = diag(S);
plot(sigmas(1:min([length(sigmas), 20])));


%%
function [IM, F] = KronFilter2(IM, levelSizes)    
    levelSize = levelSizes(1);
    IMp = zeros(floor(size(IM,1) / levelSize), floor(size(IM,1) / levelSize));
    Fp  = zeros(levelSize, levelSize);
    for i=1:size(IMp)
        for j=1:size(IMp)
            IMp(i,j) = mean(IM((i-1)*levelSize+1:i*levelSize, (j-1)*levelSize+1:j*levelSize), 'all');
            Fp = Fp + IM((i-1)*levelSize+1:i*levelSize, (j-1)*levelSize+1:j*levelSize);
        end
    end
    Fp = Fp ./ (i*j);
    if length(levelSize) > 1
        [IMs, Fs] = KronFilter(IMp, levelSizes(2:end));
        IM = cell(length(IMs) + 1, 1);
        F = cell(length(IMs) + 1, 1);
        for i=1:length(IM)
            IM{i} = IMs{i};
            F{i} = Fs{i};
        end
        IM{length(IM)} = IMp;
        F{length(F)} = Fp;
    else
        F = {Fp};
        IM = {IMp};
    end
end

