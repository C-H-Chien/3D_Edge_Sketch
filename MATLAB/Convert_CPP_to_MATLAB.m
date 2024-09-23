% fileID = fopen('pairededge6n16_T-less_FixNumofEdge.txt','r');
% formatSpec = '%f';
% A = fscanf(fileID,formatSpec);
% B = reshape(A,50,size(A,1)/50);
% B=B+1;
% paired_edge = B';
% idx_less0 = find(paired_edge<0);
% paired_cpp = paired_edge;
% paired_cpp(idx_less0) = paired_cpp(idx_less0) + 1;
% paired_edge=paired_cpp;

% load('F:\pipeline_epipolar_wedges_MATLAB\GenData\paired_tless_6n13_25.mat')
% % % % % % % % paired_edge = correct_pair;
% % % % % % % % 
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\TO_Edges_ICL-NUIM_ofkt1.mat');
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\imageArray_ICL-NUIM_ofkt1.mat')
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\R_matrix_ICL-NUIM_ofkt1.mat')
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\T_matrix_ICL-NUIM_ofkt1.mat')
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\K_ICL-NUIM_ofkt1.mat')

% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\TO_Edges_10.mat');
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\imageArray_10.mat')
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\R_matrix_10.mat')
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\T_matrix_10.mat')
% % % % % % % % load('F:\pipeline_epipolar_wedges_MATLAB\MATData\K_10.mat')

load('F:\pipeline_epipolar_wedges_MATLAB\MATData\TO_Edges_25.mat');
load('F:\pipeline_epipolar_wedges_MATLAB\MATData\imageArray_25.mat')
load('F:\pipeline_epipolar_wedges_MATLAB\MATData\R_matrix_25.mat')
load('F:\pipeline_epipolar_wedges_MATLAB\MATData\T_matrix_25.mat')
load('F:\pipeline_epipolar_wedges_MATLAB\MATData\K_25.mat')
% % % % % % % % load('bucket_rgbd_dataset_freiburg3_cabinet_validation.mat')

params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 13;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*5;
params.ICL_DATA                 = 0;
params.multiK                   = 1;

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

TO_Edges_HYPO1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
TO_Edges_HYPO2 = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};
Gamma1s = [];

for(i = 1 : size(paired_edge,1))
    edgel_HYPO1 = TO_Edges_HYPO1(paired_edge(i,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge(i,2), :);
    if(params.multiK == 1)
        pt1 = invK1 * [edgel_HYPO1(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
        pt2 = invK2 * [edgel_HYPO2(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
    else
        pt1 = invK * [edgel_HYPO1(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
        pt2 = invK * [edgel_HYPO2(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
    end
    pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
    %pts_meters = [edgel_HYPO1(1, 1), edgels_HYPO2(finalPairIndx, 1); edgel_HYPO1(1, 2), edgels_HYPO2(finalPairIndx, 2)];
    Gamma = linearTriangulation(2, pts_meters, R21, T21);
    Gamma1s = [Gamma1s, Gamma];
end

figure;
plot3(Gamma1s(1,:), -Gamma1s(2,:), -Gamma1s(3,:), 'b.','MarkerSize',6);
axis equal;
% xlim([-2 2])
% ylim([-2 2])
% zlim([-2 2])
title '3D reconstruction result'