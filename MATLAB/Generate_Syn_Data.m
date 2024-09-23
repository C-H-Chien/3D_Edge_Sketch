%> Code Description: 
%
%> (c) LEMS, Brown University
%> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
%> Last Modified: Feb 08, 2024

clear;
close all;

rng(0);
N       = 3;     %> Number of views
nPoints = 20;     %> Number of points per view
dataset_dir = 'synthcurves_dataset\spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object\';

%> Import all curve data from N views
all_R = zeros(3,3,N);
all_T = zeros(3,1,N);
nViews = round(unifrnd(1, 99, [1, N]));
pick_Img_Points = zeros(nPoints, 2, N);
pick_Scene_Points = zeros(nPoints, 3, N);
for n = 1:N
    %> read data from the synthetic curve dataset
    [img_curves, scene_curves, pose] = readCurveSyntheticDataset(nViews(n), dataset_dir, 0, 1);
    
    %> fix the 2d img points indices from the first view
    if n == 1
        imgPts_indx = round(unifrnd(1, size(img_curves.points, 1), [nPoints, 1]));
        K = pose.K;
    end
    
    %> 2D image points, 3D scene points, and absolute R, T
    collect_third_order_edges_data{n} = [img_curves.points, img_curves.tangents];
    pick_Img_Points(:,:,n) = img_curves.points(imgPts_indx, :);
    pick_Scene_Points(:,:,n) = scene_curves.points(imgPts_indx, :);
    all_R(:,:,n) = pose.R;
    all_T(:,n)   = pose.T;
end

%> Compute relative poses
R1 = all_R(:,:,1);
R2 = all_R(:,:,2);
R3 = all_R(:,:,3);
T1 = all_T(:,1);
T2 = all_T(:,2);
T3 = all_T(:,3);

%> Test
Gamma = pick_Scene_Points(10,:,1)';
gamma = K*(R1 * Gamma + T1);
gamma = gamma(1:2,:) ./ gamma(3); %> gamma now is equal to pick_Img_Points(10,:,1)





