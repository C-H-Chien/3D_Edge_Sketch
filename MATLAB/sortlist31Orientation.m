% Epipolar wedges method is used in this main code.
clear all;
close all;

% load all files needed
% load('MATData/TO_Edges_rgbd_dataset_freiburg3_cabinet_validation.mat');
% load('MATData/imageArray_rgbd_dataset_freiburg3_cabinet_validation.mat')
% load('MATData/R_matrix_rgbd_dataset_freiburg3_cabinet_validation1.mat')
% load('MATData/T_matrix_rgbd_dataset_freiburg3_cabinet_validation1.mat')
% load('MATData/K_rgbd_dataset_freiburg3_cabinet_validation.mat')

% load('MATData/TO_Edges_ICL-NUIM_ofkt1_26to50.mat');
% load('MATData/imageArray_ICL-NUIM_ofkt1_26to50.mat')
% load('MATData/R_matrix_ICL-NUIM_ofkt1_26to50.mat')
% load('MATData/T_matrix_ICL-NUIM_ofkt1_26to50.mat')
% load('MATData/K_ICL-NUIM_ofkt1.mat')

load('MATData/TO_Edges_10.mat');
load('MATData/imageArray_10.mat')
load('MATData/R_matrix_10.mat')
load('MATData/T_matrix_10.mat')
load('MATData/K_10.mat')
% load('GenData/ambi_tless_6n16.mat')

% load('MATData/TO_Edges_25.mat');
% load('MATData/imageArray_25.mat')
% load('MATData/R_matrix_25.mat')
% load('MATData/T_matrix_25.mat')
% load('MATData/K_25.mat')
% load('GenData/ambi_tless_6n13_25.mat')
% ambi_edge_prev = ambi_edge;

% load('MATData/TO_Edges_syn.mat');
% load('MATData/imageArray_syn.mat')
% load('MATData/R_matrix_syn.mat')
% load('MATData/T_matrix_syn.mat')
% load('MATData/K_syn.mat')

% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 16;
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

%> Convert from orientation to tangent vector for each edgel
TO_Edges_Tangents = cell(params.NUM_OF_IMGS, 1);
for i = 1:params.NUM_OF_IMGS
    TO_Edgels_Orientations = collect_third_order_edges{i,1}(:,3);

    Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
    for j = 1:size(TO_Edgels_Orientations, 1)
        Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
            sin(TO_Edgels_Orientations(j,1))];
    end
    TO_Edges_Tangents{i,1} = Tgts;
end

%> Array of view indices
view_idx = [];
for i = 1:params.NUM_OF_IMGS; view_idx = [view_idx; i]; end
validation_view_indices = view_idx;
validation_view_indices([params.HYPO1_VIEW_INDX],:) = [];

%> Fetch edgels of the two hypothesis views
TO_Edges_HYPO1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
TO_Edges_HYPO2 = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};
Gamma1s        = [];

hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
rows      = size(hypo_img2, 1);
cols      = size(hypo_img2, 2);

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

finalEdgePair    = [];
paired_edge      = [];
collect_third_order_edges1 = collect_third_order_edges;

for vi = 1:size(validation_view_indices, 1)
    % get the validation view
    VALID_INDX = validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    % get edges in the validation view
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    [ore_list1bar, ore_list3_sorted, epipole_pix_view1, epipole_pix_view3]= ...
        getOreList31(TO_Edges_HYPO1, TO_Edges_VALID, R_matrix, T_matrix, params, K, VALID_INDX);
    ore_list3_sortedidx = ore_list3_sorted(:,2);
    sorted31Hypo1edge = TO_Edges_VALID(ore_list3_sortedidx,:);
    collect_third_order_edges1{VALID_INDX} = sorted31Hypo1edge;

%     table1 = table(sorted31Hypo1edge);
%     mytable = table1{:,:};
%     writematrix(mytable,['F:\pipeline_epipolar_wedges_MATLAB\sortedList\Edge_',num2str(VALID_INDX),'to6.txt'],'Delimiter','tab')
end