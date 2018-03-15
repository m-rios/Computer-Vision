% clear; clc; close all;
% I1 = imread('tsukuba/scene1.row3.col1.ppm');
% I2 = imread('tsukuba/scene1.row3.col2.ppm');
% I3 = imread('tsukuba/scene1.row3.col3.ppm');
% I4 = imread('tsukuba/scene1.row3.col4.ppm');
% I5 = imread('tsukuba/scene1.row3.col5.ppm');
% 
% % figure
% % imshow(I1);
% % [xl, yl] = getpts;
% %     
% % figure
% % imshow(I5);
% % [xr, yr] = getpts;
% 
% %% Exercise 1
% load('triplets.mat');
% disparity = abs(all(:,1) - all(:,3));
% disp_map = (disparity - min(disparity))/(max(disparity)-min(disparity)); %normalise to disparity range
% 
% disp_map = mean(reshape(disp_map, 3, []));
% 
% imshow(I1);
% hold on;
% plot(all(:,1), all(:,2), '*m');
% 
% imshow(I1);
% % fill(reshape(all(:,1), [], 3), reshape(all(:,1), [], 3), 'r');
% colormap gray;
% f = fill(reshape(all(:,1),3,[]), reshape(all(:,2),3,[]),disp_map);
% 
% % distance = sqrt((clear_points(:,1) - clear_points(:,3)).^2 + (clear_points(:,2) - clear_points(:,4)).^2);
% 
% 
% %% Ex1 - Ratio comparisson with ground truth
% 
% 
% truth = imread('tsukuba/truedisp.row3.col3.pgm');
% truth_map = unique(truth);
% truth_map = cast(truth_map(2:8), 'double');
% truth_map(3) = []; %3
% truth_map(4) = []; %5
% truth_ratios = sort(truth_map / max(truth_map))'
% 
% disp_ratios = sort(mean(reshape(disparity, 3, []))/max(mean(reshape(disparity, 3, []))));
% disp_ratios(1) = []
% 
% 
% 
% %% Exercise 2
% kernel_5 = ones(5,5);
% kernel_9 = ones(9,9);
% kernel_13 = ones(13,13);
% 
% %ad
% disp_ad_5 = area_based(I1, I5, kernel_5, 'ad', 60);
% disp_ad_9 = area_based(I1, I5, kernel_9, 'ad', 60);
% disp_ad_13 = area_based(I1, I5, kernel_13, 'ad', 60);
% 
% %sd
% disp_sd_5 = area_based(I1, I5, kernel_5, 'sd', 60);
% disp_sd_9 = area_based(I1, I5, kernel_9, 'sd', 60);
% disp_sd_13 = area_based(I1, I5, kernel_13, 'sd', 60);
% 
% %ncc
% % disp_ncc_5 = area_based(I1, I5, kernel_5, 'ncc', 60);
% % disp_ncc_9 = area_based(I1, I5, kernel_9, 'ncc', 60);
% % disp_ncc_13 = area_based(I1, I5, kernel_13, 'ncc', 60);
% 
% %save images
% imwrite(disp_ad_5, 'square_ad_5.png');
% imwrite(disp_ad_5, 'square_ad_9.png');
% imwrite(disp_ad_5, 'square_ad_13.png');
% 
% imwrite(disp_sd_5, 'square_sd_5.png');
% imwrite(disp_sd_5, 'square_sd_9.png');
% imwrite(disp_sd_5, 'square_sd_13.png');
% 
% % imwrite(disp_ncc_5, 'square_ncc_5.png');
% % imwrite(disp_ncc_5, 'square_ncc_9.png');
% % imwrite(disp_ncc_5, 'square_ncc_13.png');
% 
% %compute difference with ground truth
% truth_d = cast(truth, 'double');
% norm_truth = (truth_d - min(truth_d(:)))/(max(truth_d(:) - min(truth_d(:))));
% 
% comp_ad_5 = max(max(normxcorr2(disp_ad_5, norm_truth)));
% comp_ad_9 = max(max(normxcorr2(disp_ad_9, norm_truth)));
% comp_ad_13 = max(max(normxcorr2(disp_ad_13, norm_truth)));
% 
% comp_sd_5 = max(max(normxcorr2(disp_sd_5, norm_truth)));
% comp_sd_9 = max(max(normxcorr2(disp_sd_9, norm_truth)));
% comp_sd_13 = max(max(normxcorr2(disp_sd_13, norm_truth)));
% % 
% % comp_ncc_5 = max(max(normxcorr2(disp_ncc_5, norm_truth)));
% % comp_ncc_9 = max(max(normxcorr2(disp_ncc_9, norm_truth)));
% % comp_ncc_13 = max(max(normxcorr2(disp_ncc_13, norm_truth)));
% 
% %% Exercise 3 - gaussian kernel
% sigma=1.5;
% kernel_5 = fspecial('gaussian',5, 1.5);
% kernel_9 = fspecial('gaussian',9, 3.5);
% kernel_13 = fspecial('gaussian',13, 5);
% 
% %ad
% disp_ad_5 = area_based(I1, I5, kernel_5, 'ad', 60);
% disp_ad_9 = area_based(I1, I5, kernel_9, 'ad', 60);
% disp_ad_13 = area_based(I1, I5, kernel_13, 'ad', 60);
% 
% %sd
% disp_sd_5 = area_based(I1, I5, kernel_5, 'sd', 60);
% disp_sd_9 = area_based(I1, I5, kernel_9, 'sd', 60);
% disp_sd_13 = area_based(I1, I5, kernel_13, 'sd', 60);
% 
% %ncc
% % disp_ncc_5 = area_based(I1, I5, kernel_5, 'ncc', 60);
% % disp_ncc_9 = area_based(I1, I5, kernel_9, 'ncc', 60);
% % disp_ncc_13 = area_based(I1, I5, kernel_13, 'ncc', 60);
% 
% %save images
% imwrite(disp_ad_5, 'gaussian_1_5_ad_5.png');
% imwrite(disp_ad_5, 'gaussian_1_5_ad_9.png');
% imwrite(disp_ad_5, 'gaussian_1_5_ad_13.png');
% 
% imwrite(disp_sd_5, 'gaussian_1_5_sd_5.png');
% imwrite(disp_sd_5, 'gaussian_1_5_sd_9.png');
% imwrite(disp_sd_5, 'gaussian_1_5_sd_13.png');
% 
% % imwrite(disp_ncc_5, 'gaussian_1_5_ncc_5.png');
% % imwrite(disp_ncc_5, 'gaussian_1_5_ncc_9.png');
% % imwrite(disp_ncc_5, 'gaussian_1_5_ncc_13.png');
% 
% %compute difference with ground truth
% truth_d = cast(truth, 'double');
% norm_truth = (truth_d - min(truth_d(:)))/(max(truth_d(:) - min(truth_d(:))));
% 
% comp_ad_5 = max(max(normxcorr2(disp_ad_5, norm_truth)));
% comp_ad_9 = max(max(normxcorr2(disp_ad_9, norm_truth)));
% comp_ad_13 = max(max(normxcorr2(disp_ad_13, norm_truth)));
% 
% comp_sd_5 = max(max(normxcorr2(disp_sd_5, norm_truth)));
% comp_sd_9 = max(max(normxcorr2(disp_sd_9, norm_truth)));
% comp_sd_13 = max(max(normxcorr2(disp_sd_13, norm_truth)));
% % 
% % comp_ncc_5 = max(max(normxcorr2(disp_ncc_5, norm_truth)));
% % comp_ncc_9 = max(max(normxcorr2(disp_ncc_9, norm_truth)));
% % comp_ncc_13 = max(max(normxcorr2(disp_ncc_13, norm_truth)));

%% Exercise 4
clear; clc; close all;
% [~,keys1,loc1] = sift('gray1.pgm');
% [~,keys2,loc2] = sift('gray5.pgm');

[num, or_loc1, or_loc2, matches] = match('gray1.pgm', 'gray5.pgm');


for i = 1:length(or_loc1)
    if matches(i) > 0
        loc1(i,:) = [or_loc1(i,2),or_loc1(i,1)];
        loc2(i,:) = [or_loc2(matches(i),2),or_loc2(matches(i),1)];
    end
end
loc1 = int32(ceil(loc1));
loc2 = int32(ceil(loc2));

loc1( ~any(loc1,2), : ) = [];  %rows
loc2( ~any(loc2,2), : ) = [];  %rows

%% Compare squared difference with ground truth
clear; clc; close all;

load('sift_features');

I1 = imread('tsukuba/scene1.row3.col1.ppm');
I5 = imread('tsukuba/scene1.row3.col5.ppm');

max_y_disp = sort(abs(loc1(2,:) - loc2(2,:)));
n_y_sifted = sum(abs(loc1(2,:) - loc2(2,:)) ~=0);

sift_disp = double(abs(loc1(:,1)-loc2(:,1)));
norm_sift_disp = (sift_disp - min(sift_disp))/(max(sift_disp) - min(sift_disp));

truth = imread('tsukuba/truedisp.row3.col3.pgm');
truth_d = cast(truth, 'double');

kernel_13 = fspecial('gaussian',13, 5);
kernel_13_lin = kernel_13(7,:).*ones(13,13);

%ad
disp_ad_13_gauss = area_based(I1, I5, kernel_13, 'ad', 60);
disp_ad_13_gauss_lin = area_based(I1, I5, kernel_13_lin, 'ad', 60);

norm_truth_sift = (truth_d - min(sift_disp(:)))/(max(sift_disp(:) - min(sift_disp(:))));
norm_truth_2d = (truth_d - min(disp_ad_13_gauss(:)))/(max(disp_ad_13_gauss(:) - min(disp_ad_13_gauss(:))));
norm_truth_lin = (truth_d - min(disp_ad_13_gauss_lin(:)))/(max(disp_ad_13_gauss_lin(:) - min(disp_ad_13_gauss_lin(:))));


for i=1:length(sift_disp)
    x = loc1(i,1); y = loc1(i,2); 
    comp_sift(i) = abs(norm_truth_sift(y,x)-norm_sift_disp(i));
    comp_gauss_2d(i) = abs(norm_truth_2d(y,x)-disp_ad_13_gauss(y,x));
    comp_gauss_lin(i) = abs(norm_truth_lin(y,x)-disp_ad_13_gauss_lin(y,x));   
end
