function [IM, F] = KronFilter(IM, levelSizes)    
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
    if length(levelSizes) > 1
        [IMs, Fs] = KronFilter(IMp, levelSizes(2:end));
        IM = cell(length(IMs) + 1, 1);
        F = cell(length(IMs) + 1, 1);
        for i=1:length(IMs)
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

