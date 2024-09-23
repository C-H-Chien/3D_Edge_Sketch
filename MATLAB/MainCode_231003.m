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


% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 1;
params.HYPO2_VIEW_INDX          = 2;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.001*5;
params.ICL_DATA                 = 0;
invK = inv(K);

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

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, invK);

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
for edge_idx = 1:size(TO_Edges_HYPO1, 1)
    if mod(edge_idx, 10) == 0,  fprintf(". %d",edge_idx); end
    if mod(edge_idx, 500) == 0, fprintf("\n"); end
    edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
    if edgel_HYPO1(1,1) < 10 || edgel_HYPO1(1,1) > cols-10 || ...
            edgel_HYPO1(1,2) < 10 || edgel_HYPO1(1,2) > rows-10
        continue;
    end
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
    angle_range1 = abs(rad2deg(atan(ore_list2_sorted(1,1))) - ...
        rad2deg(atan(ore_list2_sorted(size(ore_list2_sorted,1),1))));
    % get the range for single epipolar wedge in hypo 2
    range1       = params.PERCENTAGE_FOR_EPIPOLE*angle_range1;
    % get the epipolar wedge range for the edge in hypo1 on hypo2
    thresh_ore1  = tan(deg2rad(rad2deg(atan(ore_list1bar(edge_idx,1)))...
                   -range1));
    thresh_ore2  = tan(deg2rad(rad2deg(atan(ore_list1bar(edge_idx,1)))...
                   +range1));
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
            params.HYPO1_VIEW_INDX, K, VALID_INDX);
        % find epipolar angle range for validation view
        angle_range2 = abs(rad2deg(atan(ore_list31_sorted(1,1))) - ...
            rad2deg(atan(ore_list31_sorted(size(ore_list31_sorted,1),1))));
        % get the range for single epipolar wedge in vali
        range2       = params.PERCENTAGE_FOR_EPIPOLE*angle_range2;
        % get the epipolar wedge range for the edge in hypo1 on vali
        angle31_1 = rad2deg(atan(ore_31bar))-range2;
        angle31_2 = rad2deg(atan(ore_31bar))+range2;
        thresh_ore31_1  = tan(deg2rad(angle31_1));
        thresh_ore31_2  = tan(deg2rad(angle31_2));
        % get edges in vali fall inside the epipolar wedge
        idxlist_31      = find(ore_list31_sorted(:,1) >= thresh_ore31_1 & ...
            ore_list31_sorted(:,1) <= thresh_ore31_2);
        vali_idx     = ore_list31_sorted(idxlist_31,2);
        edgels_31    = TO_Edges_VALID(vali_idx, :);
        [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgel ...
             (edgel_HYPO1, edgels_HYPO2, R_matrix, T_matrix, params, VALID_INDX, K);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the orientation list between hypo2 and validation view
        [ore_32bar, ore_list32_sorted, epipole_pix_view2_3, epipole_pix_view3_2] = ...
            getOreList_Vali(edgels_HYPO2, TO_Edges_VALID, R_matrix, T_matrix, ...
            params.HYPO2_VIEW_INDX, K, VALID_INDX);
        % find epipolar angle range for validation view
        angle_range3 = abs(rad2deg(atan(ore_list32_sorted(1,1))) - ...
            rad2deg(atan(ore_list32_sorted(size(ore_list32_sorted,1),1))));
        % get the range for single epipolar wedge in vali
        range3       = params.PERCENTAGE_FOR_EPIPOLE*angle_range3;
        % investigate all the potential pairs
        % "form" the quadrilateral
        allquad = [];
        isparal = ones(size(edgels_HYPO2,1),1);
        for(hypo2idx = 1:size(edgels_HYPO2,1))
            % get the epipolar wedge range for the edge in hypo2 on vali
            angle32_1 = rad2deg(atan(ore_32bar(hypo2idx,1)))-range3;
            angle32_2 = rad2deg(atan(ore_32bar(hypo2idx,1)))+range3;
            thresh_ore32_1 = tan(deg2rad(angle32_1));
            thresh_ore32_2 = tan(deg2rad(angle32_2));
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
            anglediff    = [abs(angle31_1 - angle32_1); ...
                            abs(angle31_1 - angle32_2); ...
                            abs(angle31_2 - angle32_1); ...
                            abs(angle31_2 - angle32_2)];
            maxanglediff = max(anglediff);
            minanglediff = min(anglediff);
            if (maxanglediff <=30 && minanglediff <=30 )
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
        if max_support_val < params.MAX_NUM_OF_SUPPORT_VIEWS
            continue;
        else
            max_support_indx = find(rep_count == max_support_val);
            if size(max_support_indx, 1) > 1
                indx = unique_supported_indices(max_support_indx,1);
                finalPairIndx = getFinalHypothesisPair(Epipolar_Coeffs, edgels_HYPO2, indx, params);
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
                
                finalEdgePair = [finalEdgePair; [edgel_HYPO1(1, 1:2), edgels_HYPO2(finalPairIndx, 1:2)]];

                pt1 = invK * [edgel_HYPO1(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                pt2 = invK * [edgels_HYPO2(finalPairIndx, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
                Gamma = linearTriangulation(2, pts_meters, R21, T21);
                Gamma1s = [Gamma1s, Gamma];

                %{
                figure;
                subplot(1,2,1);
                hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
                hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
                imshow(uint8(hypo_img1)); hold on;
                plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
                plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mx', 'MarkerSize',7, 'LineWidth', 2);
                title 'hypothesis view 1'
                subplot(1,2,2);
                imshow(uint8(hypo_img2)); hold on;
                cols = size(hypo_img2, 2);
                yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
                yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
                plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
                line([1, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 1); hold on;
                plot(edgels_HYPO2(finalPairIndx, 1), edgels_HYPO2(finalPairIndx, 2), 'mx', 'MarkerSize',7, 'LineWidth', 2);
                title 'hypothesis view 2'
                %}
            end
        end
    end
    hypo2_idx   = sort(hypo2_idx);
    paired_edge = [paired_edge; edge_idx, hypo2_idx(finalPairIndx), supported_link(finalPairIndx,:)];

end

fprintf("\n> Finished!\n");
toc
figure;
plot3(Gamma1s(1,:), Gamma1s(2,:), Gamma1s(3,:), 'b.');
axis equal;
ylim([-2 2])
xlim([-2 2])
zlim([-2 2])

figure;
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
imshow(uint8(hypo_img1)); hold on;
plot(TO_Edges_HYPO1(4000:5000, 1), TO_Edges_HYPO1(4000:5000, 2), 'c.');
set(gcf,'color','w');