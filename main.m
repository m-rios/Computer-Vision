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

kernel = ones(7,7);
disp = area_based(I1,I2, kernel, 'ad');

