% Epipolar wedges method is used in this main code.
clear all; 
close all;

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
candidate_num = [];
for edge_idx = 2013:size(TO_Edges_HYPO1, 1)
    if mod(edge_idx, 10) == 0,  fprintf(". %d",edge_idx); end
    if mod(edge_idx, 500) == 0, fprintf("\n"); end
    edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
    if edgel_HYPO1(1,1) < 10 || edgel_HYPO1(1,1) > cols-10 || ...
            edgel_HYPO1(1,2) < 10 || edgel_HYPO1(1,2) > rows-10
        continue;
    end

    angle_range1 = abs(ore_list2_sorted(1,1) - ...
        ore_list2_sorted(size(ore_list2_sorted,1),1));
    range1       = params.PERCENTAGE_FOR_EPIPOLE*angle_range1;
    thresh_ore1  = ore_list1bar(edge_idx,1)-range1;
    thresh_ore2  = ore_list1bar(edge_idx,1)+range1;
    idxlist      = find(ore_list2_sorted(:,1) >= thresh_ore1 & ...
                       ore_list2_sorted(:,1) <= thresh_ore2);
    if(isempty(idxlist))
        continue;
    end
    hypo2_idx    = ore_list2_sorted(idxlist,2);
    edgels_HYPO2 = TO_Edges_HYPO2((hypo2_idx), :);
    candidate_num = [candidate_num; edge_idx, size(hypo2_idx,1)];

    % % correct edgels_HYPO2 here

    [edgels_HYPO2_corrected, edgels_HYPO1_corrected] = edgelsHYPO2correct(edgel_HYPO1, ...
        edgels_HYPO2, R21, T21, K, F, params);

    
    edgels_HYPO1_final = edgels_HYPO1_corrected;
    edgels_HYPO2_final = edgels_HYPO2_corrected;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  visualize the epipolar wedge, candidate gamma 2 and corrected gamma 2  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure(1)
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    subplot(1,2,1);
    hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
    hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
    imshow(uint8(hypo_img1));
    hold on;
    plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
    hypo31_1 = edgel_HYPO1(1,1:2)' - epipole_pix_view1(1:2,:);
    b31_1    = hypo31_1(2,:)./hypo31_1(1,:);
    c31_1    = edgel_HYPO1(1,2) - b31_1 * edgel_HYPO1(1,1);
    y31_min1 = b31_1*1 + c31_1;
    y31_max1 = b31_1*cols + c31_1;
    hold on;
    line([1, cols], [y31_min1, y31_max1], 'Color', 'g', 'LineWidth', 1);
    hold on;
%     plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'Color',"#D95319",'Marker',"x", 'MarkerSize',7, 'LineWidth', 2);
    plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mx', 'MarkerSize',7, 'LineWidth', 1);
    title 'hypothesis view 1'
    
    subplot(1,2,2);
    imshow(uint8(hypo_img2)); 
    hold on;
    cols = size(hypo_img2, 2);
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
    plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
    line([1, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 1); hold on;
    hold on;
    plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
    hypo31_1 = edgels_HYPO2(1,1:2)' - epipole_pix_view2(1:2,:);
    b31_1    = hypo31_1(2,:)./hypo31_1(1,:);
    c31_1    = edgels_HYPO2(1,2) - b31_1 * edgels_HYPO2(1,1);
    y31_min1 = b31_1*1 + c31_1;
    y31_max1 = b31_1*cols + c31_1;
    hold on;
    line([1, cols], [y31_min1, y31_max1], 'Color', 'y', 'LineWidth', 1);
    hypo31_1 = edgels_HYPO2(size(edgels_HYPO2,1),1:2)' - epipole_pix_view2(1:2,:);
    b31_1    = hypo31_1(2,:)./hypo31_1(1,:);
    c31_1    = edgels_HYPO2(size(edgels_HYPO2,1),2) - b31_1 * edgels_HYPO2(size(edgels_HYPO2,1),1);
    y31_min1 = b31_1*1 + c31_1;
    y31_max1 = b31_1*cols + c31_1;
    hold on;
    line([1, cols], [y31_min1, y31_max1], 'Color', 'y', 'LineWidth', 1);
    plot(edgels_HYPO2(:, 1), edgels_HYPO2(:, 2), 'mx', 'MarkerSize',7, 'LineWidth', 1);
    plot(edgels_HYPO2_corrected(:, 1), edgels_HYPO2_corrected(:, 2), 'gx', 'MarkerSize',7, 'LineWidth', 1);
    
    title 'hypothesis view 2'
    %}
end

fprintf("\n> Finished!\n");
figure(2);
histogram(candidate_num(:,2), 'BinWidth',1)
ylabel ('Number of edges in HYPO1', 'FontSize',20)
xlabel ('Number of edges in HYPO2 being paired with ONE edge in HYPO1', 'FontSize',20)
title ('Distribution of Number of edges in HYPO2 being paired with ONE edge in HYPO1, bin size = 1', 'FontSize',20)

table_num_of_hypo = [];
for(Hypo1 = 1:50)
    params.HYPO1_VIEW_INDX          = Hypo1;
    for(Hypo2 = 1:50)
        if(Hypo2 == Hypo1)
            continue
        end
        params.HYPO2_VIEW_INDX          = Hypo2;
        TO_Edges_HYPO1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
        TO_Edges_HYPO2 = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};
        candidate_num = [];
        [ore_list1bar, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= ...
    getOreList(TO_Edges_HYPO1, TO_Edges_HYPO2, R_matrix, T_matrix, params, K);
        for edge_idx = 1:size(TO_Edges_HYPO1, 1)
            edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
            if edgel_HYPO1(1,1) < 10 || edgel_HYPO1(1,1) > cols-10 || ...
                    edgel_HYPO1(1,2) < 10 || edgel_HYPO1(1,2) > rows-10
                continue;
            end
        
            angle_range1 = abs(ore_list2_sorted(1,1) - ...
                ore_list2_sorted(size(ore_list2_sorted,1),1));
            range1       = params.PERCENTAGE_FOR_EPIPOLE*angle_range1;
            thresh_ore1  = ore_list1bar(edge_idx,1)-range1;
            thresh_ore2  = ore_list1bar(edge_idx,1)+range1;
            idxlist      = find(ore_list2_sorted(:,1) >= thresh_ore1 & ...
                               ore_list2_sorted(:,1) <= thresh_ore2);
            if(isempty(idxlist))
                continue;
            end
            hypo2_idx    = ore_list2_sorted(idxlist,2);
            edgels_HYPO2 = TO_Edges_HYPO2(sort(hypo2_idx), :);
            candidate_num = [candidate_num; edge_idx, size(hypo2_idx,1)];
        
            % % correct edgels_HYPO2 here
        
        % %     [edgels_HYPO2_corrected, edgels_HYPO1_corrected] = edgelsHYPO2correct(edgel_HYPO1, ...
        % %         edgels_HYPO2, R21, T21, K, F, params);
        % % 
        % %     
        % %     edgels_HYPO1_final = edgels_HYPO1_corrected;
        % %     edgels_HYPO2_final = edgels_HYPO2_corrected;
        end
        table_num_of_hypo = [table_num_of_hypo; Hypo1, size(collect_third_order_edges{Hypo1},1), Hypo2, size(collect_third_order_edges{Hypo2},1), sum(candidate_num(:,2))];
    end
    
end