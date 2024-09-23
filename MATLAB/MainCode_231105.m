% Epipolar wedges method is used in this main code.
clear all; 
% close all;

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
% thinner_pair = load('GenData/wedge0006.mat').paired_edge;

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
params.ANGLE_FOR_EPIPOLE        = 0.01;
params.ICL_DATA                 = 0;
params.multiK                   = 1;
params.cols                     = size(double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX})), 2);

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
GammaTangent_3D= [];

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

% % % % load('lostidx.mat')
for edge_idx = 1:size(TO_Edges_HYPO1, 1)
% % % thinnestpair = load('GenData/verythin_test_oldcorrect.mat').paired_edge;
% % % for loss_i = 1:size(thinnestpair, 1)
% % %     edge_idx = thinnestpair(loss_i,1);
%     if (isempty(find(thinnestpair(:,1) == edge_idx)) == 0)
%         continue;
%     end
    if mod(edge_idx, 100) == 0,  fprintf(". %d",edge_idx); end
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
    range1       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range1;
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

    % % correct edgels_HYPO2 here

    [edgels_HYPO2_corrected] = edgelsHYPO2correct_1(Epipolar_Coeffs, edgels_HYPO2);    
%     [edgels_HYPO2_corrected, edgels_HYPO1_corrected] = edgelsHYPO2correct(Epipolar_Coeffs, edgels_HYPO2, edgel_HYPO1, R_matrix, T_matrix, params, K);
    edgels_HYPO2_final = edgels_HYPO2_corrected;
    edgels_HYPO1_final = edgel_HYPO1;
%     edgels_HYPO1_final = edgels_HYPO1_corrected;
%     edgels_HYPO2_final = edgels_HYPO2;
    

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
            getOreList_Vali(edgels_HYPO1_final, TO_Edges_VALID, R_matrix, T_matrix, ...
            params.HYPO1_VIEW_INDX, K, VALID_INDX, params);
        % find epipolar angle range for validation view
        angle_range2 = abs(ore_list31_sorted(1,1)- ...
            ore_list31_sorted(size(ore_list31_sorted,1),1));
        % get the range for single epipolar wedge in vali
        range2       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range2;
        % get the epipolar wedge range for the edge in hypo1 on vali
        thresh_ore31_1 = ore_31bar-range2;
        thresh_ore31_2 = ore_31bar+range2;
        idxlist_31      = find(ore_list31_sorted(:,1) >= thresh_ore31_1 & ...
            ore_list31_sorted(:,1) <= thresh_ore31_2);
        vali_idx     = ore_list31_sorted(idxlist_31,2);
        edgels_31    = TO_Edges_VALID(vali_idx, :);
        [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgel ...
             (edgels_HYPO1_final, edgels_HYPO2_final, R_matrix, T_matrix, params, VALID_INDX, K);
%         [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgelCorrected ...
%              (edgels_HYPO1_final, edgels_HYPO2_final, R_matrix, T_matrix, params, VALID_INDX, K);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the orientation list between hypo2 and validation view
        [ore_32bar, ore_list32_sorted, epipole_pix_view2_3, epipole_pix_view3_2] = ...
            getOreList_Vali(edgels_HYPO2_final, TO_Edges_VALID, R_matrix, T_matrix, ...
            params.HYPO2_VIEW_INDX, K, VALID_INDX, params);
        % find epipolar angle range for validation view
        angle_range3 = abs(ore_list32_sorted(1,1) - ...
            ore_list32_sorted(size(ore_list32_sorted,1),1));
        % get the range for single epipolar wedge in vali
        range3       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range3;
        % investigate all the potential pairs
        % "form" the quadrilateral
        allquad = [];
        isparal = ones(size(edgels_HYPO2_final,1),1);
        for(hypo2idx = 1:size(edgels_HYPO2_final,1))
            % get the epipolar wedge range for the edge in hypo2 on vali
            thresh_ore32_1 = ore_32bar(hypo2idx,1)-range3;
            thresh_ore32_2 = ore_32bar(hypo2idx,1)+range3;
            % get edges in vali fall inside the epipolar wedge
            idxlist_32 = find(ore_list32_sorted(:,1) >= thresh_ore32_1 & ...
                ore_list32_sorted(:,1) <= thresh_ore32_2);
            vali_idx1 = ore_list32_sorted(idxlist_32,2);
            edgels_32 = TO_Edges_VALID(vali_idx1, :);
            % find the edges falls inside both two wedges
            quad_idx  = intersect(vali_idx, vali_idx1);
            quad_edge = TO_Edges_VALID(quad_idx, :);
            commonedgeidx{hypo2idx} = quad_idx;
            edge32{hypo2idx} = edgels_32;
            allquad = [allquad, quad_idx'];
            anglediff    = [abs(thresh_ore31_1 - thresh_ore32_1); ...
                            abs(thresh_ore31_1 - thresh_ore32_2); ...
                            abs(thresh_ore31_2 - thresh_ore32_1); ...
                            abs(thresh_ore31_2 - thresh_ore32_2)];
            maxanglediff = max(anglediff);
            minanglediff = min(anglediff);
            if (maxanglediff <=30)
                isparal(hypo2idx,1) = 0;
            end
        end
%         if(isempty(allquad))
%             indice = [];
%             supported_indices_stack       = [supported_indices_stack; indice];
%             supported_link                = [supported_link, zeros(size(commonedgeidx,2),1)];
%         else
            Supports_orient = getSupportedEdgels_Orientation(Tangents_VALID, ...
                          reproj_edge_tgt_gamma3, params, commonedgeidx, isparal);
%         supported_edges_hypo2vali{vi} = Supports_orient.indices;
            supported_indices_stack       = [supported_indices_stack; Supports_orient.indices];
            supported_link                = [supported_link, Supports_orient.link];
%         end
    end
    if isempty(supported_indices_stack)
        continue;
    else
        %> Find number of repetitive support indices
        unique_supported_indices = unique(supported_indices_stack(:,1));
        rep_count = arrayfun(@(x)length(find(supported_indices_stack(:,1) == x)), unique_supported_indices, 'Uniform', false);
        rep_count = cell2mat(rep_count);
        [max_support_val, ~] = max(rep_count);
%{
figure(1)
len = 3;
tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
edgel_tgt1  = tgt1(edge_idx, 1:2);
hypo2_idx   = sort(hypo2_idx);
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};
subplot(1,2,1);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
imshow(uint8(hypo_img1));
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
x1 = edgel_HYPO1(1,1)+len*edgel_tgt1(1,1);
x2 = edgel_HYPO1(1,1)-len*edgel_tgt1(1,1);
y1 = edgel_HYPO1(1,2)+len*edgel_tgt1(1,2);
y2 = edgel_HYPO1(1,2)-len*edgel_tgt1(1,2);
hold on;
plot([x1 x2], [y1 y2], 'm', 'LineWidth', 1);
hold on;
hypo31_all = edgel_HYPO1(1,1:2)' - epipole_pix_view1(1:2,:);
b31_all    = hypo31_all(2,:)./hypo31_all(1,:);
c31_all    = edgel_HYPO1(1,2) - b31_all * edgel_HYPO1(1,1);
y31_min1 = b31_all*1 + c31_all;
y31_max1 = b31_all*cols + c31_all;
hold on;
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mx', 'MarkerSize',7, 'LineWidth', 1);
hold on;
plot(edgels_HYPO1_final(:, 1), edgels_HYPO1_final(:, 2), 'gx', 'MarkerSize',7, 'LineWidth', 1);
abs_R1 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C1 = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
          T_matrix(:,params.HYPO2_VIEW_INDX);
abs_R2 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
          T_matrix(:,params.HYPO1_VIEW_INDX);
R21    = abs_R2 * abs_R1';
T21    = abs_R2 * (abs_C1 - abs_C2);
% Calculation 2
% R21 =  R_matrix(:,:,new_idx    ) * R_matrix(:,:,reference_idx)';
% T21 = -R_matrix(:,:,new_idx    ) * R_matrix(:,:,reference_idx)' * T_matrix(:,:,reference_idx) + T_matrix(:,:,new_idx    );
%> Calculate Essential matrix
T_x = @(T)[0,      -T(3,1),  T(2,1); ...
           T(3,1),  0,      -T(1,1); ...
          -T(2,1),  T(1,1),  0];
E   = T_x(T21) * R21;
% Calculate fundamental matrix
if(params.multiK == 1)
    K1 = K(:,:,params.HYPO2_VIEW_INDX);
    K2 = K(:,:,params.HYPO1_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
    F12   = invK2'*E*invK1;
else
    invK = inv(K);
    F12   = invK'*E*invK;
end
for(idxH1 = 1: size (edgels_HYPO1_final,1))
    coeffs = F12 * [edgels_HYPO2_final(idxH1,1:2)'; 1];
%     a2_line = -coeffs(1,1)./coeffs(2,1);
%     b2_line = -1;
%     c2_line = -coeffs(3,1)./coeffs(2,1);
    yMin = -coeffs(3,1)./coeffs(2,1);
    yMax = (-coeffs(3,1) - coeffs(1,1)*cols) ./ coeffs(2,1);
    line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 0.5);
end
title 'hypothesis view 1'

subplot(1,2,2);
imshow(uint8(hypo_img2));
title 'hypothesis view 2'
hold on;
plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 0.5);
% y = k1*x + b1
k1 = tan(thresh_ore1/180*pi);
b1 = epipole_pix_view2(2,1) - k1*epipole_pix_view2(1,1);
yMax1 = k1*cols + b1;
% y = k2*x + b2
k2 = tan(thresh_ore2/180*pi);
b2 = epipole_pix_view2(2,1) - k2*epipole_pix_view2(1,1);
yMax2 = k2*cols + b2;
hold on;
line([0, cols], [b1, yMax1], 'Color', 'y', 'LineWidth', 0.5);
hold on;
line([0, cols], [b2, yMax2], 'Color', 'y', 'LineWidth', 0.5);
plot(edgels_HYPO2(:, 1), edgels_HYPO2(:, 2), 'mx', 'MarkerSize',7, 'LineWidth', 1);
plot(edgels_HYPO2_final(:, 1), edgels_HYPO2_final(:, 2), 'gx', 'MarkerSize',7, 'LineWidth', 1);
len = 1;
for(Indx = 1: size(edgels_HYPO2_final,1))
    edgel_tgt2  = tgt2(hypo2_idx(Indx), 1:2);
    hold on;
    x1 = edgels_HYPO2_final(Indx, 1)+len*edgel_tgt2(1,1);
    x2 = edgels_HYPO2_final(Indx, 1)-len*edgel_tgt2(1,1);
    y1 = edgels_HYPO2_final(Indx, 2)+len*edgel_tgt2(1,2);
    y2 = edgels_HYPO2_final(Indx, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'm', 'LineWidth', 0.25);
    hold on;
%     x1 = edgels_HYPO2(Indx, 1)+len*edgel_tgt2(1,1);
%     x2 = edgels_HYPO2(Indx, 1)-len*edgel_tgt2(1,1);
%     y1 = edgels_HYPO2(Indx, 2)+len*edgel_tgt2(1,2);
%     y2 = edgels_HYPO2(Indx, 2)-len*edgel_tgt2(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.5);
%     hold on;
%     plot(edgels_HYPO2_final(Indx, 1), edgels_HYPO2_final(Indx, 2), 'gx', 'MarkerSize',7, 'LineWidth', 1);
end
%}
        if max_support_val < params.MAX_NUM_OF_SUPPORT_VIEWS
            continue;
        else
            max_support_indx = find(rep_count == max_support_val);
            if size(max_support_indx, 1) > 1
                indx = unique_supported_indices(max_support_indx,1);
                finalPairIndx = getFinalHypothesisPair(Epipolar_Coeffs, edgels_HYPO2_final, indx, params);
                if(finalPairIndx ~= 0)
                    ambi_edge = [ambi_edge; edge_idx];
                end
            else
                finalPairIndx = unique_supported_indices(max_support_indx,1);
            end
            if finalPairIndx == 0
                continue;
            else
                cnt_sup = cnt_sup+1;
                sup_num{cnt_sup} = [edge_idx, edge_idx; unique_supported_indices, rep_count];

%                 finalEdgePair = [finalEdgePair; [edgel_HYPO1(1, 1:2), edgels_HYPO2_final(finalPairIndx, 1:2)]];
                finalEdgePair = [finalEdgePair; [edgels_HYPO1_final, edgels_HYPO2_final(finalPairIndx, 1:2)]];

%                 if(params.multiK == 1)
%                     pt1 = invK1 * [edgel_HYPO1(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
%                     pt2 = invK2 * [edgels_HYPO2_final(finalPairIndx, 1:2)'; 1]; pt2 = pt2(1:2, 1);
%                 else
%                     pt1 = invK * [edgel_HYPO1(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
%                     pt2 = invK * [edgels_HYPO2_final(finalPairIndx, 1:2)'; 1]; pt2 = pt2(1:2, 1);
%                 end

                if(params.multiK == 1)
                    pt1 = invK1 * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK2 * [edgels_HYPO2_final(finalPairIndx, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                else
                    pt1 = invK * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK * [edgels_HYPO2_final(finalPairIndx, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                end


                pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
                Gamma = linearTriangulation(2, pts_meters, R21, T21);
                Gamma1s = [Gamma1s, Gamma];
                GammaTgt_3D     = get3DOrientation(edgels_HYPO1_final(1, :),edgels_HYPO2_final(finalPairIndx, :), K, params, R21);
                GammaTangent_3D = [GammaTangent_3D, GammaTgt_3D];

                %{
                figure;
                %
                len = 3;
                tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
                edgel_tgt1  = tgt1(edge_idx, 1:2);
                hypo2_idx   = sort(hypo2_idx);
                tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};
                edgel_tgt2  = tgt2(hypo2_idx(finalPairIndx), 1:2);
                subplot(1,2,1);
                hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
                hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
                imshow(uint8(hypo_img1)); 
                hold on;
                plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
                hold on;
                x1 = edgel_HYPO1(1,1)+len*edgel_tgt1(1,1);
                x2 = edgel_HYPO1(1,1)-len*edgel_tgt1(1,1);
                y1 = edgel_HYPO1(1,2)+len*edgel_tgt1(1,2);
                y2 = edgel_HYPO1(1,2)-len*edgel_tgt1(1,2);
                hold on;
%                 plot([x1 x2], [y1 y2], 'Color',"#D95319", 'LineWidth', 1.5);
                hypo31_1 = edgel_HYPO1(1,1:2)' - epipole_pix_view1(1:2,:);
                b31_1    = hypo31_1(2,:)./hypo31_1(1,:);
                c31_1    = edgel_HYPO1(1,2) - b31_1 * edgel_HYPO1(1,1);
                y31_min1 = b31_1*1 + c31_1;
                y31_max1 = b31_1*cols + c31_1;
                hold on;            
%                 line([1, cols], [y31_min1, y31_max1], 'Color', 'y', 'LineWidth', 0.5);
                hold on;
                plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
                title 'hypothesis view 1'
                
                subplot(1,2,2);
                imshow(uint8(hypo_img2)); 
                title 'hypothesis view 2'
                hold on;
                plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
                len = 30;
                % for(Indx = finalPairIndx: size(edgels_HYPO2_final,1))
                Indx = finalPairIndx;
                    edgel_tgt2  = tgt2(hypo2_idx(Indx), 1:2);
%                     imshow(uint8(hypo_img2)); 
                    hold on;
                    cols = size(hypo_img2, 2);
                    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
                    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
                    
                    line([0, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 0.5); 
    %                 hold on;
    %                 x1 = edgels_HYPO2_final(Indx, 1)+len*edgel_tgt2(1,1);
    %                 x2 = edgels_HYPO2_final(Indx, 1)-len*edgel_tgt2(1,1);
    %                 y1 = edgels_HYPO2_final(Indx, 2)+len*edgel_tgt2(1,2);
    %                 y2 = edgels_HYPO2_final(Indx, 2)-len*edgel_tgt2(1,2);
    %                 hold on;
    %                 plot([x1 x2], [y1 y2], 'm', 'LineWidth', 1);
                    hold on;
                    plot(edgels_HYPO2(Indx, 1), edgels_HYPO2(Indx, 2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
                    hold on;
                    x1 = edgels_HYPO2(Indx, 1)+len*edgel_tgt2(1,1);
                    x2 = edgels_HYPO2(Indx, 1)-len*edgel_tgt2(1,1);
                    y1 = edgels_HYPO2(Indx, 2)+len*edgel_tgt2(1,2);
                    y2 = edgels_HYPO2(Indx, 2)-len*edgel_tgt2(1,2);
                    hold on;
                    plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.5);                
                    hold on;
                    plot(edgels_HYPO2_final(Indx, 1), edgels_HYPO2_final(Indx, 2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
                    hold on;
                    plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
%                 end
                title 'hypothesis view 2'
                
                %}
            end
        end
    end
    hypo2_idx   = sort(hypo2_idx);
    paired_edge = [paired_edge; edge_idx, hypo2_idx(finalPairIndx), supported_link(finalPairIndx,:)];
% BREAKPOINT
end

fprintf("\n> Finished!\n");
toc
figure(2);
plot3(Gamma1s(1,:), -Gamma1s(2,:), -Gamma1s(3,:), 'b.');
axis equal;
title '3D reconstruction result (H2 only, wedge = ± 0.01°)'
% ylim([-2 2])
% xlim([-2 2])
% zlim([-2 2])

fprintf("\n> Finished!\n");

% figure(3);
% hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
% imshow(uint8(hypo_img1)); hold on;
% plot(TO_Edges_HYPO1(1000:2000, 1), TO_Edges_HYPO1(1000:2000, 2), 'c.');
% set(gcf,'color','w');