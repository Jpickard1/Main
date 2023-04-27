%% Kronecker Dynamics
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 19, 2023

clear; close all; clc;

for i=1:70
    is = string(i);
    if i < 10; is = "0" + is; end;
    img1 = imread("2023-03-28-BJPF-s06-Stitched_t0" + is + "c1.tif");
    img2 = imread("2023-03-28-BJPF-s06-Stitched_t0" + is + "c2.tif");
    img3 = imread("2023-03-28-BJPF-s06-Stitched_t0" + is + "c3.tif");
    figure; 
    subplot(1,3,1); imagesc(img1); title(string(i));
    subplot(1,3,2); imagesc(img2);
    subplot(1,3,3); imagesc(img3);
end
