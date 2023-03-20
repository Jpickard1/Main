% Read the video
Sa = VideoReader('mouse 44_fed fast refed with motion correction for umap.mp4');
Sa_frame1 = read(Sa,[1 1508]);
Sa_frame2 = read(Sa,[1509 3017]);
Sa_frame3 = read(Sa,[3018 Inf]);

% Center the data
Data1 = double(reshape(Sa_frame1(:,:,1,:), size(Sa_frame1, 1)*size(Sa_frame1, 2), size(Sa_frame1, 4)));
nData1 = Data1 - mean(Data1, 2);

Data2 = double(reshape(Sa_frame2(:,:,1,:), size(Sa_frame2, 1)*size(Sa_frame2, 2), size(Sa_frame2, 4)));
nData2 = Data2 - mean(Data2, 2);

Data3 = double(reshape(Sa_frame3(:,:,1,:), size(Sa_frame3, 1)*size(Sa_frame3, 2), size(Sa_frame3, 4)));
nData3 = Data3 - mean(Data3, 2);

% SVD
[U1, ~, ~] = svd(nData1, 'econ');
[U2, ~, ~] = svd(nData2, 'econ');
[U3, ~, ~] = svd(nData3, 'econ');
u1 = reshape(U1(:,1), size(Sa_frame1, 1), size(Sa_frame1, 2));
u2 = reshape(U2(:,1), size(Sa_frame2, 1), size(Sa_frame2, 2));
u3 = reshape(U3(:,1), size(Sa_frame3, 1), size(Sa_frame3, 2));

% Binarize
eps = 0.015;
u1 = u1 >= eps;
u2 = u2 >= eps;
u3 = u3 >= eps;

% Plot
figure, 
imshow(u1(:,:),[])
figure, 
imshow(u2(:,:),[])
figure,
imshow(u3(:,:),[])