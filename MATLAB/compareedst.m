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
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*10;
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
validation_view_indices([params.HYPO1_VIEW_INDX; params.HYPO2_VIEW_INDX],:) = [];

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

% get the orientation list between hypo1 and hypo2
[ore_list1bar, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= ...
    getOreList(TO_Edges_HYPO1, TO_Edges_HYPO2, R_matrix, T_matrix, params, K);
tic
% pipeline starts
cnt_sup   = 0;
sup_num   = [];
ambi_edge = [];
anssss = [];
for edge_idx = 1:size(TO_Edges_HYPO1, 1)
%     if mod(edge_idx, 10) == 0,  fprintf(". %d",edge_idx); end
    if mod(edge_idx, 500) == 0, fprintf("\n"); end
    edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
    if edgel_HYPO1(1,1) < 10 || edgel_HYPO1(1,1) > cols-10 || ...
            edgel_HYPO1(1,2) < 10 || edgel_HYPO1(1,2) > rows-10
        continue;
    end
%     if(isempty(find(ambi_edge_prev == edge_idx)) == 1)
%         continue
%     end
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
%     [~, Epipolar_Coeffs, HYPO2_idx] = PairEdgeHypothesis ...
%     (F, edgel_HYPO1, TO_Edges_HYPO2);
    % find epipolar angle range for hypo 2
    angle_range1 = abs(ore_list2_sorted(1,1) - ...
        ore_list2_sorted(size(ore_list2_sorted,1),1));
    % get the range for single epipolar wedge in hypo 2
    range1       = params.PERCENTAGE_FOR_EPIPOLE*angle_range1;
    % get the epipolar wedge range for the edge in hypo1 on hypo2
    thresh_ore1  = ore_list1bar(edge_idx,1)-range1;
    thresh_ore2  = ore_list1bar(edge_idx,1)+range1;
    % get edges in hypo2 fall inside the epipolar wedge
    idxlist      = find(ore_list2_sorted(:,1) >= thresh_ore1 & ...
                       ore_list2_sorted(:,1) <= thresh_ore2);
    if(isempty(idxlist))
        continue;
    end
    hypo2_idx    = ore_list2_sorted(idxlist,2);
    edgels_HYPO2 = TO_Edges_HYPO2(sort(hypo2_idx), :);

    edgels_HYPO2_forplot = TO_Edges_HYPO2(hypo2_idx, :); 
    
    supported_indices_stack = [];
    supported_link          = [];
    
    anssss = [anssss;round((idxlist(1,1)+idxlist(size(idxlist,1),1))/2), size(TO_Edges_HYPO2,1)];
    % validate the potential pairs in all validation views
    for vi = 1:size(validation_view_indices, 1)
        % get the validation view
        VALID_INDX = validation_view_indices(vi, 1);
        VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));        
        % get edges in the validation view
        TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
        Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
        % get the orientation list between hypo1 and validation view
        [ore_31bar, ore_list31_sorted, epipole_pix_view1_3, epipole_pix_view3_1] = ...
            getOreList_Vali(edgel_HYPO1, TO_Edges_VALID, R_matrix, T_matrix, ...
            params.HYPO1_VIEW_INDX, K, VALID_INDX, params);
        % find epipolar angle range for validation view
        angle_range2 = abs(ore_list31_sorted(1,1)- ...
            ore_list31_sorted(size(ore_list31_sorted,1),1));
        % get the range for single epipolar wedge in vali
        range2       = params.PERCENTAGE_FOR_EPIPOLE*angle_range2;
        % get the epipolar wedge range for the edge in hypo1 on vali
        thresh_ore31_1 = ore_31bar-range2;
        thresh_ore31_2 = ore_31bar+range2;
        idxlist_31      = find(ore_list31_sorted(:,1) >= thresh_ore31_1 & ...
            ore_list31_sorted(:,1) <= thresh_ore31_2);
        vali_idx     = ore_list31_sorted(idxlist_31,2);
        edgels_31    = TO_Edges_VALID(vali_idx, :);
        anssss = [anssss;round((idxlist_31(1,1)+idxlist_31(size(idxlist_31,1),1))/2), size(TO_Edges_VALID,1)];
    end
end

fprintf("\n> Finished!\n");
toc
figure(1);
plot3(Gamma1s(1,:), Gamma1s(2,:), Gamma1s(3,:), 'b.');
axis equal;
% ylim([-2 2])
% xlim([-2 2])
% zlim([-2 2])

fprintf("\n> Finished!\n");

figure(2);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
imshow(uint8(hypo_img1)); hold on;
plot(TO_Edges_HYPO1(1000:2000, 1), TO_Edges_HYPO1(1000:2000, 2), 'c.');
set(gcf,'color','w');