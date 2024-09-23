% Epipolar wedges method is used in this main code.
clear all; 
% close all;

load('MATData/TO_Edges_10.mat');
load('MATData/imageArray_10.mat')
load('MATData/R_matrix_10.mat')
load('MATData/T_matrix_10.mat')
load('MATData/K_10.mat')
load('GenData/wedge0008.mat')

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
params.PERCENTAGE_FOR_EPIPOLE   = 0.005;
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

hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
rows      = size(hypo_img2, 1);
cols      = size(hypo_img2, 2);

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

finalEdgePair    = [];
% paired_edge      = [];

% get the orientation list between hypo1 and hypo2
[ore_list1bar, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= ...
    getOreList(TO_Edges_HYPO1, TO_Edges_HYPO2, R_matrix, T_matrix, params, K);
tic
% pipeline starts
cnt_sup   = 0;
sup_num   = [];
ambi_edge = [];

% load('lostidx.mat')
% for edge_idx = 1:size(TO_Edges_HYPO1, 1)
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};
figure(1);
subplot(1,2,1);
imshow(uint8(hypo_img1));
hold on;
subplot(1,2,2);
imshow(uint8(hypo_img2));
hold on;

% % % figure(2);
% % % subplot(1,2,1);
% % % imshow(uint8(hypo_img1));
% % % hold on;
% % % subplot(1,2,2);
% % % imshow(uint8(hypo_img2));
% % % hold on;
% % % 
% % % figure(3);
% % % subplot(1,2,1);
% % % imshow(uint8(hypo_img1));
% % % hold on;
% % % subplot(1,2,2);
% % % imshow(uint8(hypo_img2));
% % % hold on;
% % % 
% % % figure(4);
% % % subplot(1,2,1);
% % % imshow(uint8(hypo_img1));
% % % hold on;
% % % subplot(1,2,2);
% % % imshow(uint8(hypo_img2));
% % % hold on;
% % % 
% % % figure(5);
% % % subplot(1,2,1);
% % % imshow(uint8(hypo_img1));
% % % hold on;
% % % subplot(1,2,2);
% % % imshow(uint8(hypo_img2));
% % % hold on;
colorlist = [1 0 0;
             0 1 0;
             0 0 1;
             0 1 1;
             1 0 1;
             1 1 0;
             1 1 1;
             1 0 0;
             0 1 0];
for loss_i = 99:size(paired_edge, 1)
    fig_idx = mod(loss_i,9);
    if(fig_idx == 0)
        subplot(1,2,1);
        imshow(uint8(hypo_img1));
        hold on;
        subplot(1,2,2);
        imshow(uint8(hypo_img2));
        hold on;
    end
    color_bar = colorlist(fig_idx+1,:);%rand(1,3);
    edge_idx  = paired_edge(loss_i,1);
    hypo2_index   = paired_edge(loss_i,2);
    edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
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
len = 6;
tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
edgel_tgt1  = tgt1(edge_idx, 1:2);
hypo2_idx   = sort(hypo2_idx);
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};

% figure(fig_idx)
subplot(1,2,1);
% imshow(uint8(hypo_img1));
% hold on;
% plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 1);
hold on;
x1 = edgel_HYPO1(1,1)+len*edgel_tgt1(1,1);
x2 = edgel_HYPO1(1,1)-len*edgel_tgt1(1,1);
y1 = edgel_HYPO1(1,2)+len*edgel_tgt1(1,2);
y2 = edgel_HYPO1(1,2)-len*edgel_tgt1(1,2);
hold on;
hypo31_all = edgel_HYPO1(1,1:2)' - epipole_pix_view1(1:2,:);
b31_all    = hypo31_all(2,:)./hypo31_all(1,:);
c31_all    = edgel_HYPO1(1,2) - b31_all * edgel_HYPO1(1,1);
y31_min1 = b31_all*1 + c31_all;
y31_max1 = b31_all*cols + c31_all;
hold on;
plot([x1 x2], [y1 y2], 'Color', color_bar, 'LineWidth', 1);
hold on;
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'Color', color_bar, 'Marker', 'x','MarkerSize',8, 'LineWidth', 1);
title 'hypothesis view 1'

subplot(1,2,2);
% imshow(uint8(hypo_img2));
% title 'hypothesis view 2'
hold on;
% plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 1);
hold on;
cols = size(hypo_img2, 2);
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;

% line([0, cols], [yMin, yMax], 'Color', color_bar, 'LineWidth', 0.5);
% hold on;
plot(TO_Edges_HYPO2(hypo2_index, 1), TO_Edges_HYPO2(hypo2_index, 2), 'Color', color_bar, 'Marker', 'x','MarkerSize',8, 'LineWidth', 1);
hold on;
edgel_tgt2  = tgt2(hypo2_index, 1:2);
x1 = TO_Edges_HYPO2(hypo2_index, 1)+len*edgel_tgt2(1,1);
x2 = TO_Edges_HYPO2(hypo2_index, 1)-len*edgel_tgt2(1,1);
y1 = TO_Edges_HYPO2(hypo2_index, 2)+len*edgel_tgt2(1,2);
y2 = TO_Edges_HYPO2(hypo2_index, 2)-len*edgel_tgt2(1,2);
hold on;
plot([x1 x2], [y1 y2], 'Color', color_bar, 'LineWidth', 1);
title 'hypothesis view 2'
end

fprintf("\n> Finished!\n");
toc
% % figure(1);
% % plot3(Gamma1s(1,:), Gamma1s(2,:), Gamma1s(3,:), 'b.');
% % axis equal;
% ylim([-2 2])
% xlim([-2 2])
% zlim([-2 2])

fprintf("\n> Finished!\n");

% % figure(2);
% % hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
% % imshow(uint8(hypo_img1)); hold on;
% % plot(TO_Edges_HYPO1(1000:2000, 1), TO_Edges_HYPO1(1000:2000, 2), 'c.');
% % set(gcf,'color','w');