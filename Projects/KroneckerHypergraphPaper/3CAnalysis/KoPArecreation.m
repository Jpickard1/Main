%% KoPA Paper
%
% Auth: Joshua Pickard and Vivan
%       jpic@umich.edu
% Date: June 27, 2023
clear all; close all; clc
path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
path2data = 'C:\Joshua\MissingData\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';

%% Load Data
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
s=1; j=1;
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

D = D(1:4900,1:4900);

clearvars -except D path2data

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
itrs = 3;
maxN = 100;
vals = 2:10:maxN;
timeKSVD = zeros(length(vals), itrs);
timeBins = zeros(length(vals), itrs);
errorKSVD = zeros(length(vals),itrs);
errorBins = zeros(length(vals),itrs);
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

vals2 = vals .^2;
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
% path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
s=1; j=1;
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
if ~exist('D')
    T = readLAMMPSoutput(filePath);
    cords = getCords(T);
    D = squareform(pdist(cords));
end

clearvars -except D path2data

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

% figure;
% tiledlayout(1,4);
% nexttile;
% scree(D);
% nexttile;
% scree(Rsvd);
% nexttile;
% scree(Rnkp);
% nexttile;
% scree(Rna);

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

%% Failure Modes of this Algorithm

n = 70;
U = rand(n^2, 3); Sigma = diag(sort(rand(3,1),'descend')); V = rand(n^2, 3);
D = U * Sigma * V';

% SVD approximation
tic; [U, Sigma, V] = svds(D,1); toc
Rsvd = Sigma(1,1) * U(:,1) * V(:,1)'; % + Sigma(2,2) * U(:,2) * V(:,2)';

% NKP
tic; [B1, C1, ~] = nearestKroneckerProduct(D, [n n], [n n]); toc
Rnkp = kron(B1, C1);

% New Alg.
tic; [Dp, Fp] = KronFilter(D, n); toc
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

D = D / sum(D, 'all');
Rsvd = Rsvd / sum(Rsvd, 'all');
Rnkp = Rnkp / sum(Rnkp, 'all');
Rna = Rna / sum(Rna, 'all');

norm(D)
norm(D - Rsvd)
norm(D - Rnkp)
norm(D - Rna)


%% Higher Rank examples
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
% path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
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


%% Higher rank approximation with KronFilter

% Make data
n = 16;
M = zeros(n^2);
M(n^2/2,n^2/2) = 1;
R = bwdist(M);
A = R >= 45;
A = double(A); A = A * 10;
% A = A(1:289,1:289); n = 17;
itrs = 5;
figure; tiledlayout(itrs, 4);
% nexttile; imagesc(A); title('Data'); disp(norm(A));
Ac = cell(1,itrs);
Bc = cell(1,itrs);
Cc = cell(1,itrs);
Rc = cell(1,itrs);
for i=1:itrs
    [B, C] = KronFilter2(A, n); R = kron(B, C);
    nexttile; imagesc(A); title('Remaining Data');
    nexttile; imagesc(R); title('Approximation');
    nexttile; imagesc(B); title('B');
    nexttile; imagesc(C); title('C');
    fprintf("Data Norm: %f \t Approximation Norm: %f \t B Norm: %f \t C Norm: %f\n", norm(A), norm(R), norm(B), norm(C));
    Ac{i} = A;
    Bc{i} = B;
    Cc{i} = C;
    Rc{i} = R;
    A = A - R;
end

Rt = zeros(size(A));
for i=1:length(Rc)
    Rt = Rt + Rc{i};
end
figure; imagesc(Rt)

fs = 10;
figure; imagesc(conv2(Rt, ones(fs, fs)))
figure; imagesc(conv2(Rt, Bc{1}))



%%
B1 = rand(3,3);
C1 = rand(3,3);
B2 = rand(3,3);
C2 = rand(3,3);
A = kron(B1, C1) + kron(B2, C2);
figure;
tiledlayout(1,5);
nexttile; imagesc(A); title('Data')
nexttile; imagesc(B1); title('B1')
nexttile; imagesc(C1); title('C1')
nexttile; imagesc(B2); title('B2')
nexttile; imagesc(C2); title('C2')

figure; tiledlayout(itrs, 4);
% nexttile; imagesc(A); title('Data'); disp(norm(A));
for i=1:itrs
    [B, C] = KronFilter2(A, 3); R = kron(B, C);
    nexttile; imagesc(A); title('Remaining Data');
    nexttile; imagesc(R); title('Approximation');
    nexttile; imagesc(B); title('B');
    nexttile; imagesc(C); title('C');
    fprintf("Data Norm: %f \t Approximation Norm: %f \t B Norm: %f \t C Norm: %f\n", norm(A), norm(R), norm(B), norm(C));
    A = A - R;
end


%% Plot of number of parameters used
n = 10080;
x = [];
p = [];
for i=1:n
    if mod(n,i) == 0
        x = [x i];
        p = [p (i^2+(n/i)^2)/(n^2)];
    end
end

figure; scatter(x,p,'.');
title('Resolution v. Parameters');
xlabel("Factor Matrix Size");
ylabel("Parameters");
set(gca, 'YScale', 'log', 'XScale', 'log');

%% Higher Rank examples
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 3;
s=1; j=1;
path2sims = path2data + SIMPARMS(s) + "\conf.";
filePath = path2sims + string(j) + ".DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
D = squareform(pdist(cords));

clearvars -except D path2data

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
Rnkp = kron(B1,C1) + kron(B2,C2) + kron(B3,C3); % + kron(B4,C4) + kron(B5,C5) + kron(B6,C6);

%% New Alg.
Ac = cell(1,itrs);
Bc = cell(1,itrs);
Cc = cell(1,itrs);
Rc = cell(1,itrs);
A = D; itrs = 2; n = 70;
Rc = cell(1,itrs);
figure; tiledlayout(itrs, 4);
An = zeros(itrs,1);
for i=1:itrs
    % A = A / norm(A);
    [B, C] = KronFilter2(A, n); R = kron(B, C);
    nexttile([1 2]); imagesc([A R]); title('Remaining Data + Estimate');
    % nexttile; imagesc(R); title('Approximation');
    nexttile; imagesc(B); title('B');
    nexttile; imagesc(C); title('C');
    An(i) = norm(A, 1);
    fprintf("Data Norm: %f \t Approximation Norm: %f \t B Norm: %f \t C Norm: %f\n", An(i), norm(R, 1), norm(B, 1), norm(C, 1));
    Ac{i} = A;
    Bc{i} = B;
    Cc{i} = C;
    Rc{i} = R;
    A = A - R; A(A<0) = 0;
end

Rna = zeros(size(A));
for i=1:length(Rc)
    Rna = Rna + (1/i) * Rc{i};
end

figure; imagesc(Rna)
figure; imagesc(Rnkp)


%%
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

%% lsqnonneg for kronecker approximation
%   A = B kron C
clear
A = sym('A_%d', [4 4])
B = sym('B_%d', [2 2])
B_vec = B(:)
% S2 = shuffleMatrix(2,2)
% S2 * A

% Compute the permutation matrix
permutation = reshape(1:16, [4 4]).';
permutation = [1 2 5 6;
               3 4 7 8;
               9 10 13 14;
               11 12 15 16];

Ap = kron2vecPerm(A)

% Rearrange the elements of A using the permutation matrix
% A_reordered = A(permutation);

lsqnonneg(rand(4,1), rand(4,1))


%%
clear
n = 3;
B = rand(n,n)
C = rand(n,n)
A = kron(B,C)

% B = sym("B_%d", [2 2]);
% C = sym("C_%d", [2 2]);
% A = kron(B,C)

Bvec = B(:);

% Create the kronecker product matrix for B
kronB = kron(eye(n*n), Bvec);

Ap = kron2vecPerm(A);
Avec = Ap(:);

% Solve for C using lsqnonneg
C_vec = lsqnonneg(kronB, Avec);

% Reshape C to the appropriate size
Cnew = reshape(C_vec, [n n])'
C




%% Reshape A and B to vectors
S = shuffleMatrix(2,2);
Ap = A * S;

Ap(:)

C

%% Tensor based Kronecker binning
clear
itrs = 3;
maxN = 10;
vals = 2:2:maxN;
timeTKPSVD = zeros(length(vals), itrs);
timeBins = zeros(length(vals), itrs);
errorTKPSVD = zeros(length(vals),itrs);
errorBins = zeros(length(vals),itrs);
for i=1:length(vals)
    n = vals(i);
    for j=1:itrs
        A = rand(n^2, n^2, n^2);
        % tic; [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);
        tic; [Bk, sigmas] = tkpsvd(A, n*ones(1,6));
        timeTKPSVD(i,j) = toc;
        tic; [Bb, Cb] = HyperKronFilter(A, n);
        timeBins(i,j) = toc;
    
        disp("===========");
        nA = norm(A, 'fro');
        disp(nA);
        disp(norm(sigmas(1) * superkron(Bk{1,1},Bk{2,1}), 'fro'));
        disp(norm(2 * superkron(Bb,Cb), 'fro'));
    
        errorTKPSVD(i,j) = norm(A - sigmas(1) * superkron(Bk{1,1},Bk{2,1}), 'fro') / nA;
        errorBins(i,j) = norm(A - 2 * superkron(Bb,Cb), 'fro') / nA;
    end
end

vals2 = vals .^2;
vals2p = repmat(vals2, [1 itrs]);
figure;
subplot(1,2,1); hold on;
scatter(vals2p, reshape(timeTKPSVD, [1 numel(timeTKPSVD)]), 'x', 'r');
scatter(vals2p, reshape(timeBins, [1 numel(timeBins)]), '.', 'b');
title('Run Time'); xlabel('Tensor Dimension'); ylabel('Seconds');
subplot(1,2,2); hold on;
scatter(vals2p, reshape(errorTKPSVD, [1 numel(errorTKPSVD)]), 'x', 'r');
scatter(vals2p, reshape(errorBins, [1 numel(errorBins)]), '.', 'b');
title('Reconstruction Error'); xlabel('Tensor Dimension'); ylabel('2 Norm');
legend(["KSVD", "New Alg."]);
sgtitle("Run Time and Error Analysis of Rank 1 Kronecker Approximation (3-way Tensor)");

%% Rank 1 KSVD vs SVD Approximations

maxItrs = 300;
n2 = 144; n = round(sqrt(n2));
Ac = cell(maxItrs,1);
errors = zeros(maxItrs,2);

for i=1:45
    disp(i);
    % A = rand(n2,n2);
    A = D{i,s}; A = A(1:4900,1:4900); n = 70;
    [U, Sigma, V] = svds(A, 1);
    R1 = U * Sigma * V';
    [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);
    R2 = kron(Bk, Ck);
    errors(i,1) = norm(A - R1) / norm(A);
    errors(i,2) = norm(A - R2) / norm(A);
    % Ac{i} = A;
end

% hold on;
figure; 
scatter(errors(1:i-1,1), errors(1:i-1,2), '.');

%% Write all Hip-Hop distance matrices to file
SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 15;
path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
path2data = 'C:\Joshua\MissingData\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\';
D = cell(nsamps,3);
for s=1:3
    path2sims = path2data + SIMPARMS(s) + "\conf.";
    for j=1:nsamps
        filePath = path2sims + string(j) + ".DNA";
        T = readLAMMPSoutput(filePath);
        cords = getCords(T);
        D{j,s} = squareform(pdist(cords));
        disp(j);
    end
end

save("HipHopD3_45.mat",D)
save("HipHopD3_45.mat",'D','-v7.3')

%% KSVD vs New Alg

maxItrs = 100;
n = 20;
n2 = n^2;
Ac = cell(maxItrs,1);
errors = zeros(maxItrs,2);

for i=1:maxItrs
    disp(i);
    A = rand(n2,n2);
    [Bb, Cb] = KronFilter(A, n);
    [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);
    R1 = kron(Bb{1}, Cb{1});
    R2 = kron(Bk, Ck);
    errors(i,1) = norm(A - R1) / norm(A);
    errors(i,2) = norm(A - R2) / norm(A);
end

figure; 
scatter(errors(1:i-1,1), errors(1:i-1,2), '.');


%% KSVD vs New Alg. on Permuted Matrix
clear; clc; close all;
maxItrs = 1000;
n = 8;
n2 = n^2;
A = rand(n2,n2);
A = triu(A) - diag(diag(A)); A = A + A';
Ac = cell(maxItrs,1);
errors = zeros(maxItrs,2);

for i=1:maxItrs
    disp(i);
    Ac{i} = A;
    [Bb, Cb] = KronFilter(A, n);
    [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);

    R1 = kron(Bb{1}, Cb{1});
    R2 = kron(Bk, Ck);
    
    errors(i,1) = norm(A - R1) / norm(A);
    errors(i,2) = norm(A - R2) / norm(A);

    % Swap a row and a column of A
    if rand() > 0.5
        r1 = randi([1 n2]);
        r2 = randi([1 n2]);
        A([r1 r2],:) = A([r2 r1], :);
    else
        c1 = randi([1 n2]);
        c2 = randi([1 n2]);
        A([c1 c2],:) = A([c2 c1], :);
    end
end

figure; 
plot(errors(:,1), errors(:,2));
xlabel("KronFilter Error");
ylabel("NKP Error");

% scatter(errors(1:i-1,1), errors(1:i-1,2), '.');

%% KSVD vs SVD Rank 1 on Permuted Matrix
clear; clc; close all;
maxItrs = 1000;
n = 5;
n2 = n^2;
A = rand(n2,1) * rand(1, n2);
A = kron(rand(n), rand(n))
% A = rand(n2,n2);
% A = triu(A) - diag(diag(A)); A = A + A';
Ac = cell(maxItrs,1);
errors = zeros(maxItrs,2);

for i=1:maxItrs
    disp(i);
    Ac{i} = A;
    % [Bb, Cb] = KronFilter(A, n);
    [U, Sigmas, V] = svds(A,1);
    [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);

    R1 = U * Sigmas * V';
    R2 = kron(Bk, Ck);
    
    errors(i,1) = norm(A - R1) / norm(A);
    errors(i,2) = norm(A - R2) / norm(A);

    % Swap a row and a column of A
    if rand() > 0.5
        r1 = randi([1 n2]);
        r2 = randi([1 n2]);
        A([r1 r2],:) = A([r2 r1], :);
    else
        c1 = randi([1 n2]);
        c2 = randi([1 n2]);
        A([c1 c2],:) = A([c2 c1], :);
    end
end

figure; 
plot(errors(:,1), errors(:,2));
xlabel("SVD Error");
ylabel("NKP Error");

figure;
subplot(2,3,1); imagesc(Ac{1}); title("Data");
subplot(2,3,2); imagesc(R1); title("SVD");
subplot(2,3,3); imagesc(R2); title("KNP");
subplot(2,3,5); imagesc(A-R1); title("SVD Error");
subplot(2,3,6); imagesc(A-R2); title("KNP Error");

norm(A)
norm(R1)
norm(R2)

norm(A - R1) / norm(A)
norm(A - R2) / norm(A)



% scatter(errors(1:i-1,1), errors(1:i-1,2), '.');
%% Verifying nearestKroneckerProduct.m
%   result: It seems good at what it does
clear

m = 4;
n = 6;
A = sym('A_%d%d',[n m])
% n1 = round(sqrt(n)); n2 = n1;
n1 = 3; n2 = 2;
m1 = round(sqrt(m)); m2 = m1;

R = reshape( permute( reshape( A, [ m2, m1, n2, n1 ] ), [ 2 4 1 3 ] ), m1 * n1, m2 * n2 );

reverse_perm = reshape(permute(reshape(R, [n2, m2, n1, m1]), [3, 1, 4, 2]), m1, n1, m2, n2);

% Reconstruct the original matrix A
Ar = reshape(reverse_perm, m1 * n1, m2 * n2);


B = rand(n);
C = rand(n);
A = kron(B,C);
[Bg, Cg] = nearestKroneckerProduct(A, [n n], [n n])

B - Bg

C - Cg

Ag = kron(Bg, Cg)
A - Ag % Don't forget to loop to see that the scaling factor of this matrix is very small 1e-15 when n = 6 for like 5 tries in a row (thats basically correct)


%% KSVD vs SVD a few comparisons
%   rank 1 matrix
%   kronecker rank 1 matrix
%   random symmetric matrix
%   rank 1 + kronecker rank 1
%   random matrix
% clear; clc; close all;
maxItrs = 10000;
n = 10;
n2 = n^2;
matrixTypes = ["Rank 1", "Kronecker Rank 1", "Rank 1 + Kronecker Rank 1", "Random Symmetric", "Random Matrix"];
figure;
tiledlayout(2,3);
for type=1:5
    if type == 1
        A = rand(n2,1) * rand(1, n2);
    elseif type == 2
        A = kron(rand(n), rand(n));
    elseif type == 3
        A = rand(n2, 1) * rand(1,n2) + kron(rand(n), rand(n));
    elseif type == 4
        A = rand(n2,n2);
        A = triu(A); A = A + A' - diag(diag(A));
    elseif type == 5
        A = rand(n2, n2);
    else
        error('Whats going on here?');
    end
    errors = zeros(maxItrs,2);    
    for i=1:maxItrs
        if mod(i, 100) == 0; disp(i); end
        [U, Sigmas, V] = svds(A,1);
        [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);
    
        R1 = U * Sigmas * V';
        R2 = kron(Bk, Ck);
        
        errors(i,1) = norm(A - R1) / norm(A);
        errors(i,2) = norm(A - R2) / norm(A);
    
        % Swap a row and a column of A
        r1 = randi([1 n2]);
        r2 = randi([1 n2]);
        if rand() > 0.5
            A([r1 r2],:) = A([r2 r1], :);
        else
            A(:,[r1 r2]) = A(:,[r2 r1]);
        end
    end
    nexttile;
    scatter(errors(:,1), errors(:,2), '.');
    xlabel("SVD Error");
    ylabel("NKP Error");
    title(matrixTypes(type));
end

%% SVD vs KSVD for Hi-C
% clear all; close all; clc
% load("HipHopD3_45.mat")

n = 70;
for s=1:3
    for i=1:15
        disp(i);
        A = D{i,s};
        A = A(1:4900,1:4900);
        
        [U, Sigmas, V] = svds(A,1);
        [Bk,Ck, ~] = nearestKroneckerProduct(A, [n n], [n n]);
        
        R1 = U * Sigmas * V';
        R2 = kron(Bk, Ck);
        
        errors((s-1)*15+i,1) = norm(A - R1) / norm(A);
        errors((s-1)*15+i,2) = norm(A - R2) / norm(A);
        
    end
end

x = 0.2:0.01:0.45;

figure;
scatter(errors(:,1), errors(:,2), '.'); hold on;
plot(x,x);
xlabel("SVD Error");
ylabel("NKP Error");
title("Simulated Hi-C (HipHop LAMMPS)")
% title(matrixTypes(type));

figure;
imagesc(D{6,1})


%  11, 3
%  13, 1
%   6, 1

%% KSVD vs Kronecker 
clear
itrs = 100;
n = 5; n2 = n^2;
errors = zeros(itrs,2);
for i=1:itrs
    disp(i);
    % Make tensor and unfolding
    A = rand(n2,n2,n2);
    Am = block3unfold(A, n);
    % Kronecker approximate
    [Bg, Cg] = KronFilter5(Am, n, n2);
    [Bk,Ck, ~] = nearestKroneckerProduct(Am, [n n2], [n n2]);
    R1 = kron(Bg, Cg);
    R2 = kron(Bk, Ck);
    % Calculate errors
    E1 = Am - R1;
    E2 = Am - R2;

    errors(i,1) = sum(abs(E1), 'all');
    errors(i,2) = sum(abs(E2), 'all');
end

figure;
plot(errors(:,1)); hold on;
plot(errors(:,2)); hold on;

%%
clear
itrs = 3;
maxN = 5;
vals = 2:1:12;
timeKSVD = zeros(length(vals), itrs);
timeBins = zeros(length(vals), itrs);
errorKSVD = zeros(length(vals),itrs);
errorBins = zeros(length(vals),itrs);
for i=1:length(vals)
    n = vals(i); n2 = n^2;
    for j=1:itrs
        disp(j);
        % Make tensor and unfolding
        A = rand(n2,n2,n2);
        Am = block3unfold(A, n);

        % Kronecker approximate
        tic;
        [Bg, Cg] = KronFilter5(Am, n, n2);
        timeBins(i,j) = toc;
        tic;
        [Bk,Ck, ~] = nearestKroneckerProduct(Am, [n n2], [n n2]);
        timeKSVD(i,j) = toc;
        
        R1 = kron(Bg, Cg);
        R2 = kron(Bk, Ck);
    
        E1 = Am - R1;
        E2 = Am - R2;
    
        errorBins(i,j) = sum(abs(E1), 'all');
        errorKSVD(i,j) = sum(abs(E2), 'all');
    end
end

vals2 = vals .^3;
vals2p = repmat(vals2, [1 itrs]);
figure;
subplot(1,2,1); hold on;
scatter(vals2p, reshape(timeKSVD, [1 numel(timeKSVD)]), 'x', 'r');
scatter(vals2p, reshape(timeBins, [1 numel(timeBins)]), '.', 'b');
title('Run Time'); xlabel('Tensor Size'); ylabel('Seconds');
subplot(1,2,2); hold on;
scatter(vals2p, reshape(errorKSVD, [1 numel(errorKSVD)]), 'x', 'r');
scatter(vals2p, reshape(errorBins, [1 numel(errorBins)]), '.', 'b');
title('Reconstruction Error'); xlabel('Matrix Dimension'); ylabel('2 Norm');
legend(["KSVD", "New Alg."]);
sgtitle("Run Time and Error Analysis of Rank 1 Kronecker Approximation");
