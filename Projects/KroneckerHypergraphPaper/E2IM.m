function IM = E2IM(E, n)
    IM = zeros(n, size(E,1));
    for e=1:size(E,1)
        IM(E(e,1),e) = 1;
        IM(E(e,2),e) = 1;
    end
    IM = unique(IM','rows')';
end

