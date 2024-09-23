%> Code Description: Given point correspondences from three views, two 2D
%                    tangents from the first two views, and the relative 
%                    poses, find the 2D tangent from the third view.
%
%> (c) LEMS, Brown University
%> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
%> Last Modified: Jan. 10th, 2023

clear;
close all;

rng(0);
N       = 10;        %> Number of views
nPoints = 1;     %> Number of points per view
dataset_dir = 'synthcurves_dataset\spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object\';

%> Import all curve data from N views
all_R = zeros(3,3,N);
all_T = zeros(3,1,N);
% nViews =[81, 90, 13, 91, 63, 11, 28, 55, 95, 96];
nViews = round(unifrnd(1, 99, [1, N]));
pickImgPts = zeros(nPoints, 2, N);
pickImgTgt = zeros(nPoints, 2, N);
pickSceTgt = zeros(nPoints, 3, N);
for n = 1:N
    %> read data from the synthetic curve dataset
    [img_curves, scene_curves, pose] = readCurveSyntheticDataset(nViews(n), dataset_dir, 0, 1);
    
    %> fix the 2d img points indices from the first view
    if n == 1
        imgPts_indx = round(unifrnd(1, size(img_curves.points, 1), [nPoints, 1]));
        K = pose.K;
    end
    
    %> 2d image points, attached tangents, and absolute R, T
    pickImgPts(:,:,n) = img_curves.points(imgPts_indx, :);
    pickImgTgt(:,:,n) = img_curves.tangents(imgPts_indx, :);
    all_R(:,:,n) = pose.R;
    all_T(:,n)   = pose.T;

    collect_third_order_edges{n} = [img_curves.points, img_curves.tangents];
    R_matrix(:,:,n) = pose.R;
    T_matrix(:,:,n) = pose.T;
    % save the img:

end

%> Compute relative poses
R1 = all_R(:,:,1);
R2 = all_R(:,:,2);
R3 = all_R(:,:,3);
T1 = all_T(:,1,1);
T2 = all_T(:,1,2);
T3 = all_T(:,1,3);
R12 = R2 * R1';
T12 = -R2 * R1' * T1 + T2;
R13 = R3 * R1';
T13 = -R3 * R1' * T1 + T3;

%> Points and tangents in meters in view 1
gamma1           = inv(K) * [pickImgPts(:,:,1), 1]';
pt_tgt_to_pixels = pickImgPts(:,:,1) + pickImgTgt(:,:,1);
pt_tgt_to_meters = inv(K) * [pt_tgt_to_pixels, 1]';
tgt1_meters      = pt_tgt_to_meters - gamma1;

%> Points and tangents in meters in view 2
gamma2           = inv(K) * [pickImgPts(:,:,2), 1]';
pt_tgt_to_pixels = pickImgPts(:,:,2) + pickImgTgt(:,:,2);
pt_tgt_to_meters = inv(K) * [pt_tgt_to_pixels, 1]';
tgt2_meters      = pt_tgt_to_meters - gamma2;

%> Find the 3D tangent in view 1 coordinate
P1   = gamma1;
P2   = gamma2;
t1   = tgt1_meters;
t2   = tgt2_meters;
n1   = cross(t1, P1);
n2   = R12' * cross(t2, P2);   
T_v1 = cross(n1,n2) ./ norm(cross(n1,n2));

%> Transform the 3D tangent from view 1 coordinate to view 3 coordinate
T_v3 = R13 * T_v1;

%> Project from 3D tangent to 2D tangent
e3     = [0;0;1];
gamma3 = inv(K) * [pickImgPts(:,:,3), 1]';
t_v3   = T_v3 - (e3'*T_v3)*gamma3;
t_v3   = t_v3 ./ norm(t_v3);

%> The resultant 2D tangent on the third image which is the same as
%  pickImgTgt(:,:,3).
t3 = t_v3(1:2, 1);



R_matrix = all_R;
T_matrix = all_T;