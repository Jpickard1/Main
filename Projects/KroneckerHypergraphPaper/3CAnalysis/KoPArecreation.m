%% KoPA Paper
%
% Auth: Joshua Pickard and Vivan
%       jpic@umich.edu
% Date: June 27, 2023

%% Load Data
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
s=1; j=1;
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

clearvars -except D

D = D(1:4900,1:4900);

%% Compute SVD
tic;
[U, Sigma, V] = svd(D);
toc1 = toc;

R1 = Sigma(1,1) * U(:,1) * V(:,1)'; % + Sigma(2,2) * U(:,2) * V(:,2)';

sigmas = diag(Sigma);
% figure; scatter(1:length(sigmas), sigmas, '.');

%% Compute KSVD
tic;
[ B1, C1, diff1 ] = nearestKroneckerProduct(D(1:4900,1:4900), [70 70], [70 70]);
toc2 = toc;
[ B2, C2, diff2 ] = nearestKroneckerProduct(diff1, [70 70], [70 70]);
[ B3, C3, diff3 ] = nearestKroneckerProduct(diff2, [70 70], [70 70]);
[ B4, C4, diff4 ] = nearestKroneckerProduct(diff3, [70 70], [70 70]);
[ B5, C5, diff5 ] = nearestKroneckerProduct(diff4, [70 70], [70 70]);
[ B6, C6, diff6 ] = nearestKroneckerProduct(diff5, [70 70], [70 70]);

K3 = kron(B1,C1) + kron(B2,C2) + kron(B3,C3) + kron(B4,C4) + kron(B5,C5) + kron(B6,C6);

%%
R2 = kron(B,C);
figure;
tiledlayout(1,3);
nexttile;
imagesc(R2);
nexttile;
imagesc(kron(B,C));
nexttile;
imagesc(kron(C,B));

%% Make an image
figure;
tiledlayout(1,3);
nexttile;
imagesc(D); title('Data');
nexttile;
imagesc(R1); title('SVD Rank 1');
nexttile;
imagesc(R2); title('KSVD Rank 1');

%% rank 6 approximation with svd

R6 = U(:,1:6) * Sigma(1:6,1:6) * V(:,1:6)';

figure;
tiledlayout(1,4);
nexttile;
imagesc(D); title('Data');
nexttile;
imagesc(R6); title('SVD Rank 6');
nexttile;
imagesc(K3); title('KSVD Rank 6');
nexttile;
imagesc(R2); title('KSVD Rank 1');


norm(D(1:4900,1:4900) - R6(1:4900,1:4900))

norm(D(1:4900,1:4900) - K3)


%% NKP vs Joshua's Binning Method



tic;
[Dp, Fp] = KronFilter(D, 70);
toc3 = toc;
R3 = kron(Dp{1}, Fp{1});


figure;
tiledlayout(1,3);
nexttile;
imagesc(D); title('Data');
nexttile;
imagesc(R3); title('Bins');
nexttile;
imagesc(R2); title('KSVD')

norm(D(1:4900,1:4900) - R2)

norm(D(1:4970,1:4970) - R3)

R3 = max(max(D)) / (max(max(R3))) * R3;

%% Plots for paper with randomized data
itrs = 15;
maxN = 40;
vals = 2:1:maxN;
timeKSVD = zeros(length(vals), j);
timeBins = zeros(length(vals), j);
errorKSVD = zeros(length(vals),j);
errorBins = zeros(length(vals),j);
for i=1:length(vals)
    n = vals(i);
    for j=1:itrs
        A = rand(n^2, n^2);
        tic; [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);
        timeKSVD(i,j) = toc;
        tic; [Bb, Cb] = KronFilter(A, n);
        timeBins(i,j) = toc;
    
        disp("===========");
        nA = norm(A);
        disp(nA);
        disp(norm(kron(Bk,Ck)));
        disp(norm(2 * kron(Bb{1},Cb{1})));
    
        errorKSVD(i,j) = norm(A - kron(Bk,Ck)          , 2) / nA;
        errorBins(i,j) = norm(A - 2 * kron(Bb{1},Cb{1}), 2) / nA;
    end
end

vals2 = vals .^2; itrs=15;
vals2p = repmat(vals2, [1 itrs]);
figure;
subplot(1,2,1); hold on;
scatter(vals2p, reshape(timeKSVD, [1 numel(timeKSVD)]), 'x', 'r');
scatter(vals2p, reshape(timeBins, [1 numel(timeBins)]), '.', 'b');
title('Run Time'); xlabel('Matrix Dimension'); ylabel('Seconds');
subplot(1,2,2); hold on;
scatter(vals2p, reshape(errorKSVD, [1 numel(errorKSVD)]), 'x', 'r');
scatter(vals2p, reshape(errorBins, [1 numel(errorBins)]), '.', 'b');
title('Reconstruction Error'); xlabel('Matrix Dimension'); ylabel('2 Norm');
legend(["KSVD", "New Alg."]);
sgtitle("Run Time and Error Analysis of Rank 1 Kronecker Approximation");

%% Rank 1 examples
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
s=1; j=1;
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

clearvars -except D

D = D(1:4900,1:4900);

% SVD approximation
tic; [U, Sigma, V] = svds(D,1); toc
Rsvd = Sigma(1,1) * U(:,1) * V(:,1)'; % + Sigma(2,2) * U(:,2) * V(:,2)';

% NKP
tic; [B1, C1, ~] = nearestKroneckerProduct(D, [70 70], [70 70]); toc
Rnkp = kron(B1, C1);

% New Alg.
tic; [Dp, Fp] = KronFilter(D, 70); toc
Rna = kron(Dp{1}, Fp{1});

figure;
tiledlayout(1,4);
nexttile;
imagesc(D); title('Original Data');
nexttile;
imagesc(Rsvd); title('SVD Estimation Rank 1');
nexttile;
imagesc(Rnkp); title('Nearest Kronecker Product');
nexttile;
imagesc(Rna); title('New Algorithm');

figure;
tiledlayout(1,4);
nexttile;
scree(D);
nexttile;
scree(Rsvd);
nexttile;
scree(Rnkp);
nexttile;
scree(Rna);

D = D / sum(D, 'all');
Rsvd = Rsvd / sum(Rsvd, 'all');
Rnkp = Rnkp / sum(Rnkp, 'all');
Rna = Rna / sum(Rna, 'all');

norm(D)
norm(D - Rsvd)
norm(D - Rnkp)
norm(D - Rna)

A = imread("testImg1.jpg");
A = sum(A,3);
D = A(1:289,1:289);

% SVD approximation
[U, Sigma, V] = svd(D);
Rsvd = Sigma(1,1) * U(:,1) * V(:,1)'; % + Sigma(2,2) * U(:,2) * V(:,2)';

% NKP
[B1, C1, ~] = nearestKroneckerProduct(D, [17 17], [17 17]);
Rnkp = kron(B1, C1);

% New Alg.
[Dp, Fp] = KronFilter(D, 17);
Rna = kron(Dp{1}, Fp{1});

figure;
tiledlayout(1,4);
nexttile;
imagesc(D); title('Original Data');
nexttile;
imagesc(Rsvd); title('SVD Estimation Rank 1');
nexttile;
imagesc(Rnkp); title('Nearest Kronecker Product');
nexttile;
imagesc(Rna); title('New Algorithm');

%% Higher Rank examples
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
s=1; j=1;
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

clearvars -except D

D = D(1:4900,1:4900);
D = D / sum(D, 'all');

% SVD approximation
[U, Sigma, V] = svds(D,6);
Rsvd = U(:,1:6) * Sigma(1:6,1:6) * V(:,1:6)';

% NKP
[ B1, C1, diff1 ] = nearestKroneckerProduct(D, [70 70], [70 70]);
[ B2, C2, diff2 ] = nearestKroneckerProduct(diff1, [70 70], [70 70]);
[ B3, C3, diff3 ] = nearestKroneckerProduct(diff2, [70 70], [70 70]);
[ B4, C4, diff4 ] = nearestKroneckerProduct(diff3, [70 70], [70 70]);
[ B5, C5, diff5 ] = nearestKroneckerProduct(diff4, [70 70], [70 70]);
[ B6, C6, ~ ] = nearestKroneckerProduct(diff5, [70 70], [70 70]);
Rnkp = kron(B1,C1) + kron(B2,C2) + kron(B3,C3) + kron(B4,C4) + kron(B5,C5) + kron(B6,C6);

% New Alg.
% disp(sum(D,'all'));
[Dp1, Fp1] = KronFilter(D    , 70);
R1 = kron(Dp1{1},Fp1{1});
R1 = (sum(D,'all')/sum(R1,'all')) * R1;
% disp(sum(R1,'all'));
diff1 = D - R1;
% disp(sum(diff1,'all'));
[Dp2, Fp2] = KronFilter(diff1, 70);
R2 = kron(Dp2{1},Fp2{1});
R2 = (sum(diff1,'all')/sum(R2,'all')) * R2;
diff2 = diff1 - R2;
[Dp3, Fp3] = KronFilter(diff2, 70);
R3 = kron(Dp3{1},Fp3{1});
R3 = (sum(diff2,'all')/sum(R3,'all')) * R3;
diff3 = diff2 - R3;
[Dp4, Fp4] = KronFilter(diff3, 70);
R4 = kron(Dp4{1},Fp4{1});
R4 = (sum(diff3,'all')/sum(R4,'all')) * R4;
diff4 = diff3 - R4;
[Dp5, Fp5] = KronFilter(diff4, 70);
R5 = kron(Dp5{1},Fp5{1});
R5 = (sum(diff4,'all')/sum(R5,'all')) * R5;
diff5 = diff4 - R5;
[Dp6, Fp6] = KronFilter(diff5, 70);
R6 = kron(Dp6{1},Fp6{1});
R6 = (sum(diff4,'all')/sum(R6,'all')) * R6;
diff6 = diff5 - R6;
Rna = R1 + R2 + R3 + R4 + R5 + R6;

R1 = kron(Dp1{1},Fp1{1});
R2 = kron(Dp2{1},Fp2{1});
R3 = kron(Dp3{1},Fp3{1});
R4 = kron(Dp4{1},Fp4{1});
R5 = kron(Dp5{1},Fp5{1});
R6 = kron(Dp6{1},Fp6{1});

figure;
tiledlayout(1,4);
nexttile; imagesc(D); title('Data');
nexttile; imagesc(R1); title('Estimation 1');
nexttile; imagesc(D-R1); title('Error 1');
nexttile; imagesc(R2); title('Estimation 2');


sum(R1, 'all')


figure;
tiledlayout(1,4);
nexttile;
imagesc(D); title('Original Data');
nexttile;
imagesc(Rsvd); title('SVD Estimation (Rank 6)');
nexttile;
imagesc(Rnkp); title('Nearest Kronecker Product (Rank 6)');
nexttile;
imagesc(Rna); title('New Algorithm');

%%
figure; tiledlayout(1,7);
nexttile; imagesc(D);
nexttile; imagesc(diff1);
nexttile; imagesc(diff2);
nexttile; imagesc(diff3);
nexttile; imagesc(diff4);
nexttile; imagesc(diff5);
nexttile; imagesc(diff6);

%% Recursive strategy
% New Alg.
figure;
subplot(2,7,1); imagesc(D); title('Original Data');
[Dp, Fp] = KronFilter(D    , 70); diff1 = abs(D     - kron(Dp{1}, Fp{1})); Rna = kron(Dp{1}, Fp{1});
subplot(2,7,2); imagesc(diff1); title('Data 2');
subplot(2,7,8); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 1');
[Dp, Fp] = KronFilter(diff1, 70); diff2 = abs(diff1 - kron(Dp{1}, Fp{1})); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,7,3); imagesc(diff1); title('Data 3');
subplot(2,7,9); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 2');
[Dp, Fp] = KronFilter(diff2, 70); diff3 = abs(diff2 - kron(Dp{1}, Fp{1})); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,7,4); imagesc(diff1); title('Data 4');
subplot(2,7,10); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 3');
[Dp, Fp] = KronFilter(diff3, 70); diff4 = abs(diff3 - kron(Dp{1}, Fp{1})); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,7,5); imagesc(diff1); title('Data 5');
subplot(2,7,11); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 4');
[Dp, Fp] = KronFilter(diff4, 70); diff5 = abs(diff4 - kron(Dp{1}, Fp{1})); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,7,6); imagesc(diff1); title('Data 6');
subplot(2,7,12); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 5');
[Dp, Fp] = KronFilter(diff5, 70); diff6 = abs(diff5 - kron(Dp{1}, Fp{1})); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,7,7); imagesc(diff1); title('Data 7');
subplot(2,7,13); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 6');

figure; tiledlayout(1,7);
nexttile; imagesc(D);
nexttile; imagesc(diff1);
nexttile; imagesc(diff2);
nexttile; imagesc(diff3);
nexttile; imagesc(diff4);
nexttile; imagesc(diff5);
nexttile; imagesc(diff6);

A = imread("testImg1.jpg");
A = sum(A,3);
D = A(1:289,1:289); D = D / sum(D, 'all');
figure;
subplot(2,4,1); imagesc(D); title('Original Data');
[Dp, Fp] = KronFilter(D    , 17);
R = kron(Dp{1}, Fp{1}); R = R / sum(R, 'all');
diff1 = abs(D - R); Rna = kron(Dp{1}, Fp{1});
subplot(2,4,2); imagesc(diff1); title('Data 2');
subplot(2,4,5); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 1');
[Dp, Fp] = KronFilter(diff1, 17);
R = kron(Dp{1}, Fp{1}); R = R / sum(R, 'all');
diff2 = abs(D - R); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,4,3); imagesc(diff2); title('Data 3');
subplot(2,4,6); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 2');
[Dp, Fp] = KronFilter(diff2, 17);
R = kron(Dp{1}, Fp{1}); R = R / sum(R, 'all');
diff3 = abs(D - R); Rna = Rna + kron(Dp{1}, Fp{1});
subplot(2,4,4); imagesc(diff3); title('Data 4');
subplot(2,4,7); imagesc(kron(Dp{1}, Fp{1})); title('Estimation 3');
subplot(2,4,8); imagesc(Rna); title('Rank 4 Estimation');


figure;
[Dp, Fp] = KronFilter(D    , 17);
E = kron(Dp{1}, Fp{1});
diff1 = abs(D - E); Rna = kron(Dp{1}, Fp{1});
subplot(1,3,1); imagesc(D); title('Data');
subplot(1,3,2); imagesc(E); title('Estimation');
subplot(1,3,3); imagesc(diff1); title('Error');
