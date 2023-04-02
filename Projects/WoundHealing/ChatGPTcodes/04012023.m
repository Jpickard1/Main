% This code was initially written by ChatGPT and then modified by Joshua
%
% PROMPT: Please write MATLAB code to do cell segmentation based on
%         thresholds
%
% Auth: Joshua Pickard and ChatGPT
%       jpic@umich.edu
% Date: April 1, 2023

%% Joshua's Code
clear; close all; clc;
v = VideoReader('2023-03-24-triangle-stitched.avi');
e = 10;
E = cell(v.NumFrames, 1);
for fi=1:v.NumFrames
    f = read(v,fi);
    img = f(100:size(f,1)-100,50:size(f,1)-50,:);
    % Convert the image to grayscale
    img_gray = rgb2gray(img);
    
    % Kronecker Product
    % img_gray = kron(double(img_gray), ones(e, e));

    % Edge detection
    E{fi} = edge(img_gray);
    disp(fi);
end

%% 

figure;
for fi=1:length(E)
    f = read(v,fi);
    img = f(100:size(f,1)-100,50:size(f,1)-50,:);
    imagesc(img);
    hold on;
    [x, y] = find(E{fi});
    scatter(x, y);
    hold off;
    pause(0.10)
end



%% Identify edges present in most cells
C = zeros(size(img_gray));
for fi=1:length(E)
    [x, y] = find(E{fi} > 0);
    for i=1:length(x)
        C(x(i), y(i)) = C(x(i), y(i)) + 1;
    end
end
Ct = (C > 50);
figure; imagesc(Ct * 100)




%% ChatGPT Code
% Load the microscope image
% img = imread('microscope_image.tif');

img = f(100:size(f,1)-100,50:size(f,1)-50,:); figure; imagesc(img);

% Convert the image to grayscale
img_gray = rgb2gray(img);


% Kronecker Product
e = 10;
img_gray = kron(double(img_gray), ones(e, e));

E = edge(img_gray);
figure; imagesc(E);

% Apply a Gaussian filter to smooth the image and reduce noise
%for k=2:10
%    img_smooth = imgaussfilt(img_gray, k);
%    figure; imagesc(img_smooth); title("Kernel: " + string(k));
%end

img_smooth = imgaussfilt(img_gray, 3);

% Set the threshold value
thresh_value = mean(mean(img_smooth));

% Segment the image using thresholding
img_thresh = img_smooth > thresh_value;

figure; imagesc(img_thresh);

% Remove small objects (noise) from the image
img_thresh = bwareaopen(img_thresh, 50);

% Display the segmented image
imshow(img_thresh);
title('Segmented Cells');

% Calculate the number of cells
cc = bwconncomp(img_thresh);
num_cells = cc.NumObjects;
disp(['Number of cells detected: ', num2str(num_cells)]);

%% Removing stictched noise effects

% Load the stitched image
% img = imread('stitched_image.tif');

% Convert the image to grayscale
img_gray = rgb2gray(img); figure; imagesc(img_gray);

% Apply a median filter to remove salt and pepper noise
img_median = medfilt2(img_gray, [3 3]); figure; imagesc(img_median);

% Threshold the image to create a binary mask
% mask = img_median > graythresh(img_median); figure; imagesc(mask);
mask = img_median > mean(mean(img_median)); figure; imagesc(mask);

% Remove objects smaller than a certain size
mask = bwareaopen(mask, 50);

% Fill small holes in the mask
mask = imfill(mask, 'holes');

% Apply a morphological closing operation to remove small gaps in the objects
se = strel('disk', 2);
mask_closed = imclose(mask, se);

% Remove objects larger than a certain size
mask_final = bwareafilt(mask_closed, [0, 1000]);

% Apply the mask to the original image to remove the noise
img_clean = bsxfun(@times, img, cast(mask_final,class(img)));

% Display the cleaned image
imshow(img_clean);

