%% April 12, 2023
%
%   IMAGE PROCESSING PIPELINE:
%       1. Detect image boundary of stitches
%       2. For each frame in the the stitched image
%           3. Preproccess individual frames
%           4. Count cells within the frames
%           5. Add results into constructed image
%
% Auth: Joshua B. Pickard
%       jpic@umich.edu
% Date: April 12, 2023

%% ChatGPT: write cell segmentation code for a .tif image in matlab (channel 1)
clear all; close all; clc;

figure; c = 5; % Number of columns in figure / steps in process
% Read in the image
img = imread('2023-03-28-BJPF-s06-Stitched_t001c3.tif');
subplot(1,c,1); imagesc(img); title('Initial Image');

% Convert the image to grayscale
gray_img = rgb2gray(img);
subplot(1,c,2); imagesc(gray_img); title('Gray Scale');

sigma = 3;
gray_img = imgaussfilt(gray_img, sigma);
subplot(1,c,3); imagesc(gray_img); title('Gaussian Filter');
% gray_img = medfilt2(gray_img, [3 3]);
% subplot(1,c,3); imagesc(gray_img); title('Median Filter');

% Threshold the image
bw_img = zeros(size(gray_img));
for col=1:3
    cidxs = round((col-1) * (size(gray_img, 2)/3) + 1):round(col * (size(gray_img, 2)/3));
    for r = 1:7
        ridxs = round((r-1) * (size(gray_img, 1)/7) + 1):round((r * size(gray_img, 1)/7));
        gray_sub_img = gray_img(ridxs, cidxs);
        threshold = 1.7 * graythresh(gray_sub_img);
        bw_sub_img = imbinarize(gray_sub_img,threshold);
        bw_img(ridxs, cidxs) = double(bw_sub_img);

        edge_img = edge(gray_sub_img, 'Canny');

        % Combine the binary image and edge image
        combined_img = bw_sub_img & edge_img;
        
        % Label the connected components in the combined image
        [label_img, num_labels] = bwlabel(combined_img);
        
        % Display the results
        figure;
        subplot(1,2,1);
        imshow(gray_sub_img);
        title('Original Image');
        subplot(1,2,2);
        imshow(label2rgb(label_img));
        title('Connected Components');
    end
end
subplot(1,c,4); imagesc(bw_img); title('Binary Image');

%%
% Remove large objects
% Compute the complement of the binary image
complement_img = imcomplement(bw_img);

% Remove small connected components using bwareaopen
% Apply bwareafilt to retain small connected components
filtered_img = bwareafilt(bw_img, [0, 500]);
subplot(1,c,5); imagesc(filtered_img); title('Small Components Only');

% Identify Islands
% Label the image and get the region properties
labeled_img = bwlabel(bw_img);
props = regionprops(labeled_img,'Area','BoundingBox','Centroid');

% Create a figure to display the segmented cells
figure;
imshow(img);
hold on;

% Loop through the regions and plot them on the image
for i = 1:length(props)
    bbox = props(i).BoundingBox;
    x = bbox(1);
    y = bbox(2);
    w = bbox(3);
    h = bbox(4);
    rectangle('Position',[x,y,w,h],'EdgeColor','r','LineWidth',2);
end

%%

%% ChatGPT: write cell segmentation code for a .tif image in matlab (channel 1)
clear all; close all; clc;

figure; c = 5; % Number of columns in figure / steps in process
% Read in the image
img = imread('2023-03-28-BJPF-s06-Stitched_t001c1.tif');
subplot(1,c,1); imagesc(img); title('Initial Image');

% Convert the image to grayscale
gray_img = rgb2gray(img);
subplot(1,c,2); imagesc(gray_img); title('Gray Scale');

gray_img = medfilt2(gray_img, [4 4]);
subplot(1,c,3); imagesc(gray_img); title('Median Filter');

% Threshold the image
threshold = graythresh(gray_img);
bw_img = imbinarize(gray_img,threshold);
subplot(1,c,4); imagesc(bw_img); title('Binary Image');

% Remove large objects
% Compute the complement of the binary image
complement_img = imcomplement(bw_img);

% Remove small connected components using bwareaopen
% Apply bwareafilt to retain small connected components
filtered_img = bwareafilt(bw_img, [0, 500]);
subplot(1,c,5); imagesc(filtered_img); title('Small Components Only');

% Identify Islands
% Label the image and get the region properties
labeled_img = bwlabel(bw_img);
props = regionprops(labeled_img,'Area','BoundingBox','Centroid');

% Create a figure to display the segmented cells
figure;
imshow(img);
hold on;

% Loop through the regions and plot them on the image
for i = 1:length(props)
    bbox = props(i).BoundingBox;
    x = bbox(1);
    y = bbox(2);
    w = bbox(3);
    h = bbox(4);
    rectangle('Position',[x,y,w,h],'EdgeColor','r','LineWidth',2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trial Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


threshold = 100; % adjust threshold as needed
filtered_img = bwareaopen(complement_img, threshold);
subplot(1,c,4); imagesc(filtered_img); title('Binary Image');

% Remove small objects
bw_img = bwareaopen(bw_img, 50);

% Perform morphological operations to clean up the image
se = strel('disk',2);
bw_img = imclose(bw_img,se);
bw_img = imfill(bw_img,'holes');

figure; imagesc(bw_img);

% Label the image and get the region properties
labeled_img = bwlabel(bw_img);
props = regionprops(labeled_img,'Area','BoundingBox','Centroid');

% Create a figure to display the segmented cells
figure;
imshow(img);
hold on;

% Loop through the regions and plot them on the image
for i = 1:length(props)
    bbox = props(i).BoundingBox;
    x = bbox(1);
    y = bbox(2);
    w = bbox(3);
    h = bbox(4);
    rectangle('Position',[x,y,w,h],'EdgeColor','r','LineWidth',2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trial Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Fourier transform to analyze frequency content
f = fft2(gray_img);
f = fftshift(f);
mag = abs(f);

% Threshold the magnitude to identify seams
threshold = max(mag(:)) * 0.2; % adjust threshold as needed
seams = mag > threshold;

% Display the results
figure;
subplot(1,2,1);
imshow(img);
title('Original Image');
subplot(1,2,2);
imshow(seams);
title('Seams');


% Apply a median filter to remove the noise
% medImg = medfilt2(img);

sigma = 3;
gaussImg = imgaussfilt(img, sigma);

figure;
subplot(1,3,1); imagesc(img); title('Image');
% subplot(1,3,1); imagesc(medImg); title('Mediang Filter');
subplot(1,3,3); imagesc(gaussImg); title('Gaussian Smoothing, sigma=3');

% Convert the image to grayscale
gray_img = rgb2gray(img);



figure;
imshow(gray_img);

% Threshold the image
threshold = graythresh(gray_img);
bw_img = imbinarize(gray_img,threshold);

figure; imagesc(bw_img);

% Remove small objects
bw_img = bwareaopen(bw_img, 50);

figure; imagesc(bw_img);

% Perform morphological operations to clean up the image
se = strel('disk',2);
bw_img = imclose(bw_img,se);
bw_img = imfill(bw_img,'holes');

figure; imagesc(bw_img);

% Label the image and get the region properties
labeled_img = bwlabel(bw_img);
props = regionprops(labeled_img,'Area','BoundingBox','Centroid');

% Create a figure to display the segmented cells
figure;
imshow(img);
hold on;

% Loop through the regions and plot them on the image
for i = 1:length(props)
    bbox = props(i).BoundingBox;
    x = bbox(1);
    y = bbox(2);
    w = bbox(3);
    h = bbox(4);
    rectangle('Position',[x,y,w,h],'EdgeColor','r','LineWidth',2);
end


imgFile = "2023-03-28-BJPF-s06-Stitched_t001c1.tif";
t = Tiff(imgFile,"r");
readIn = read(t);

figure; imagesc(readIn)
C1 = readIn(:,:,1);
C2 = readIn(:,:,2);
C3 = readIn(:,:,3);


mean(C2)

figure;
subplot(1,3,1); imagesc(C1);
subplot(1,3,2); imagesc(C2);
subplot(1,3,3); imagesc(C3);


