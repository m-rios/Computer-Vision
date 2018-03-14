clear; clc; close all;
I1 = imread('tsukuba/scene1.row3.col1.ppm');
I2 = imread('tsukuba/scene1.row3.col2.ppm');
I3 = imread('tsukuba/scene1.row3.col3.ppm');
I4 = imread('tsukuba/scene1.row3.col4.ppm');
I5 = imread('tsukuba/scene1.row3.col5.ppm');

% figure
% imshow(I1);
% [xl, yl] = getpts;
%     
% figure
% imshow(I5);
% [xr, yr] = getpts;

%% Exercise 1
load('triplets.mat');
disparity = abs(all(:,1) - all(:,3));
disp_map = (disparity - min(disparity))/(max(disparity)-min(disparity)); %normalise to disparity range

disp_map = mean(reshape(disp_map, 3, []));

imshow(I1);
hold on;
plot(all(:,1), all(:,2), '*m');

imshow(I1);
% fill(reshape(all(:,1), [], 3), reshape(all(:,1), [], 3), 'r');
colormap gray;
f = fill(reshape(all(:,1),3,[]), reshape(all(:,2),3,[]),disp_map);

% distance = sqrt((clear_points(:,1) - clear_points(:,3)).^2 + (clear_points(:,2) - clear_points(:,4)).^2);


%% Ex1 - Ratio comparisson with ground truth


truth = imread('tsukuba/truedisp.row3.col3.pgm');
truth_map = unique(truth);
truth_map = cast(truth_map(2:8), 'double');
truth_map(3) = []; %3
truth_map(4) = []; %5
truth_ratios = sort(truth_map / max(truth_map))'

disp_ratios = sort(mean(reshape(disparity, 3, []))/max(mean(reshape(disparity, 3, []))));
disp_ratios(1) = []



%% Exercise 2
kernel_5 = ones(5,5);
kernel_9 = ones(9,9);
kernel_13 = ones(13,13);

%ad
disp_ad_5 = area_based(I1, I5, kernel_5, 'ad', 60);
disp_ad_9 = area_based(I1, I5, kernel_9, 'ad', 60);
disp_ad_13 = area_based(I1, I5, kernel_13, 'ad', 60);

%sd
disp_sd_5 = area_based(I1, I5, kernel_5, 'sd', 60);
disp_sd_9 = area_based(I1, I5, kernel_9, 'sd', 60);
disp_sd_13 = area_based(I1, I5, kernel_13, 'sd', 60);

%ncc
% disp_ncc_5 = area_based(I1, I5, kernel_5, 'ncc', 60);
% disp_ncc_9 = area_based(I1, I5, kernel_9, 'ncc', 60);
% disp_ncc_13 = area_based(I1, I5, kernel_13, 'ncc', 60);

%save images
imwrite(disp_ad_5, 'square_ad_5.png');
imwrite(disp_ad_5, 'square_ad_9.png');
imwrite(disp_ad_5, 'square_ad_13.png');

imwrite(disp_sd_5, 'square_sd_5.png');
imwrite(disp_sd_5, 'square_sd_9.png');
imwrite(disp_sd_5, 'square_sd_13.png');

% imwrite(disp_ncc_5, 'square_ncc_5.png');
% imwrite(disp_ncc_5, 'square_ncc_9.png');
% imwrite(disp_ncc_5, 'square_ncc_13.png');

%compute difference with ground truth
truth_d = cast(truth, 'double');
norm_truth = (truth_d - min(truth_d(:)))/(max(truth_d(:) - min(truth_d(:))));

comp_ad_5 = max(max(normxcorr2(disp_ad_5, norm_truth)));
comp_ad_9 = max(max(normxcorr2(disp_ad_9, norm_truth)));
comp_ad_13 = max(max(normxcorr2(disp_ad_13, norm_truth)));

comp_sd_5 = max(max(normxcorr2(disp_sd_5, norm_truth)));
comp_sd_9 = max(max(normxcorr2(disp_sd_9, norm_truth)));
comp_sd_13 = max(max(normxcorr2(disp_sd_13, norm_truth)));
% 
% comp_ncc_5 = max(max(normxcorr2(disp_ncc_5, norm_truth)));
% comp_ncc_9 = max(max(normxcorr2(disp_ncc_9, norm_truth)));
% comp_ncc_13 = max(max(normxcorr2(disp_ncc_13, norm_truth)));

corr = normxcorr2(norm_disp, norm_truth);
comparison = max(corr(:));

%% Exercise 3


