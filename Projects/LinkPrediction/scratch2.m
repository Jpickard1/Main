H = HAT.load('DBLP');
A = H * H';
A = full(A);
A = (A>0);
A = real(A);
n = size(A,1);
A = A - tril(A);
[vi, vj] = find(A == 1);
E = [vi vj];                    % Get edge set

pKnown = 0.8;

% Set list of link prediction functions
F = {@adamicAdar_index, @commonNeighbors_index, @hubDepressed_index, @hubPromoted_index, ...
    @jaccard_index, @leichtHolmeNewman_index, @salton_index, @sorensen_index }; %, ...
    % @averageCommuteTime_index, @randomWalkWithRestart_index};

% 2. Remove edges
E = E(randperm(size(E, 1)), :); % Randomly permute edge set

Ek = E(1:round(size(E,1) * pKnown), :);
Eu = E(round(size(E,1) * pKnown) + 1:end, :);

% Get edges to predict
Ec = nchoosek(1:n, 2);         % list all possible edges
Ec = setdiff(Ec, E, 'rows');   % remove edges that already exist
Ec = Ec(randperm(size(Ec, 1)), :);
Ep = [Ec(1:size(Eu), :); Eu];

% Set batch sizes
B = 0.2:0.2:1;
B = size(Eu,1) .* B;
% B = [20 50 100 150 200 300 400 500 size(Eu, 1)];

% round(exp(1:1:log(size(Eu, 1))))
% B(end + 1) = size(Eu, 1);

T = table();
for bi=1:length(B)
    b = B(bi);
    disp(b);
    bT = cell(length(F), 1);
    for fi=1:length(F)
        f = F{fi};
        disp(f);
        [Ei, EiNat] = bpredictNat(n, Ek, size(Eu, 1), b, f, Ep);
        acc = (size(Eu, 1) - size(setdiff(Eu, Ei, 'rows'), 1)) / size(Eu, 1); % / size(Eu, 1);
        T{fi, bi} = acc;
        nat = [];
        for i=1:length(EiNat)
            rNat = (size(Eu, 1) - size(setdiff(Eu, EiNat{i}, 'rows'), 1)) / size(EiNat{i}, 1);
            nat = [nat rNat];
        end
        NatNat{fi, bi} = nat;
    end
end
disp(T)