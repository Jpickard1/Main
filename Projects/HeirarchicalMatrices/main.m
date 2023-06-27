F = [1    0;
     0    1];
F = [1  1   1   0;
     1  0   1   1;
     1  1   0   1;
     1  1   1   1]
L = ones(2,2);
% L = [0.9 0.9; 0.9 0.9];

nplots = 4;
figure;
tiledlayout(1, nplots)

nexttile;
imagesc(F);

nexttile;
F2 = kron(eye(2,2), F);
imagesc(F2);

nexttile;
F3 = kron(eye(2,2), F2);
imagesc(F3)

nexttile;
F4 = kron(eye(2,2), F3);
imagesc(F4)

