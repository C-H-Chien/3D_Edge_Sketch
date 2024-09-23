close all
clear all

load('TO_Edges_ICL-NUIM_ofkt1_26to50.mat');
load('imageArray_ICL-NUIM_ofkt1_26to50.mat')
load('R_matrix_ICL-NUIM_ofkt1_26to50.mat')
load('T_matrix_ICL-NUIM_ofkt1_26to50.mat')
load('K_ICL-NUIM_ofkt1.mat')
load('wrong_paire_45n42_25_40.mat')
pair = wrong_pair40;
% load('correct_paire_45n42_25.mat')



% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 45-25;
params.HYPO2_VIEW_INDX          = 42-25;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.9995;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.001*5;
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

% for (idx_pair = 1:size(wrong_pair,1))
%     edgel_HYPO1 = TO_Edges_HYPO1(wrong_pair(idx_pair,1), :);
%     edgel_HYPO2 = TO_Edges_HYPO2(wrong_pair(idx_pair,2), :);
for (idx_pair = 232:size(pair,1))
    edgel_HYPO1 = TO_Edges_HYPO1(pair(idx_pair,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(pair(idx_pair,2), :);

    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;

    figure(1);
    subplot(1,2,1);
    hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
    hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
    imshow(uint8(hypo_img1));
    hold on;
    plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
    plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mo', 'MarkerSize',10, 'LineWidth', 2);
    plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mx', 'MarkerSize',10, 'LineWidth', 2);
    title 'hypothesis view 1'
    subplot(1,2,2);
    imshow(uint8(hypo_img2));
    hold on;
    cols = size(hypo_img2, 2);
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
    plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
    line([1, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 1); hold on;
    plot(edgel_HYPO2(1, 1), edgel_HYPO2(1, 2), 'mo', 'MarkerSize',10, 'LineWidth', 2);
    plot(edgel_HYPO2(1, 1), edgel_HYPO2(1, 2), 'mx', 'MarkerSize',10, 'LineWidth', 2);
    title 'hypothesis view 1'

    figure(2)
    for(i = 1:23)
        subplot(5,5,i)
        imshow(imageArray{validation_view_indices(i)});
        hold on;
        plot(collect_third_order_edges{validation_view_indices(i),1}(:, 1), ...
            collect_third_order_edges{validation_view_indices(i),1}(:, 2), 'c.', 'MarkerSize',0.2, 'LineWidth', 2);
%         if(wrong_pair(idx_pair,i+27) >0)
%             support_edge = collect_third_order_edges{validation_view_indices(i),1}(wrong_pair(idx_pair,i+27), :);
%             plot(support_edge(1, 1), support_edge(1, 2), 'r*');
%             plot(support_edge(1, 1), support_edge(1, 2), 'ro');
%             title(['support vali ', num2str(i)])
        if(pair(idx_pair,i+2) >0)
            support_edge = collect_third_order_edges{validation_view_indices(i),1}(pair(idx_pair,i+2), :);
            plot(support_edge(1, 1), support_edge(1, 2), 'rx', 'MarkerSize',10, 'LineWidth', 2);
            plot(support_edge(1, 1), support_edge(1, 2), 'ro', 'MarkerSize',10, 'LineWidth', 2);
            title(['support vali ', num2str(i)])
        else
            title(['vali ', num2str(i)])
        end
    end
end