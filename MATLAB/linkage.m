% close all
clear all
load('MATData/TO_Edges_10.mat');
load('MATData/imageArray_10.mat')
% R_matrixbeforeBA = load('MATData/R_matrix_10.mat').R_matrix;
% T_matrixbeforeBA = load('MATData/T_matrix_10.mat').T_matrix;
R_matrixafterBA  = load('MATData/R_matrix.mat').R_matrix;
T_matrixafterBA  = load('MATData/T_matrix.mat').T_matrix;
load('MATData/K_10.mat')
% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 12;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 15;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.ICL_DATA                 = 0;
params.multiK                   = 1;
params.syn                      = 0;
params.cols                     = size(double(rgb2gray(imageArray{1})), 2);
params.rows                     = size(double(rgb2gray(imageArray{1})), 1);
params.delta                    = 0.3;

hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
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


load('GenData\gamma1_ore15_delta03_6n12.mat')

load('GenData\pair_ore15_delta03_6n12.mat')

visual_idx = 1;

figure
plot3(Gamma1s(1,:), -Gamma1s(2,:), -Gamma1s(3,:), 'b.');
hold on
plot3(Gamma1s(1,visual_idx), -Gamma1s(2,visual_idx), -Gamma1s(3,visual_idx), 'rx', 'MarkerSize',10, 'LineWidth',2);
axis equal;
title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels, orientation different ≤ ', num2str(params.SUPPORT_OREN_THRESH),'°)'])


v_indices = 1:50;
v_indices(:,[params.HYPO1_VIEW_INDX; params.HYPO2_VIEW_INDX]) = [];
v_indices = [params.HYPO1_VIEW_INDX, params.HYPO2_VIEW_INDX, v_indices];

paired_edge_order(:,v_indices) = paired_edge;

figure(2)
imgidx = 1;
imgnum = 1;
while(imgnum<=6)
    if(paired_edge_order(visual_idx,imgidx) == 0)
        imgidx = imgidx+1;
        continue;
    end
    subplot(2,3,imgnum)
    TO_Edges = collect_third_order_edges{imgidx,1};
    imshow(imageArray{imgidx});
    hold on;
    plot(TO_Edges(:, 1), TO_Edges(:, 2), 'c.', 'MarkerSize',5);
    hold on;
    plot(TO_Edges(paired_edge_order(visual_idx,imgidx), 1), TO_Edges(paired_edge_order(visual_idx,imgidx), 2), 'rx', 'MarkerSize',10, 'LineWidth',2);
    title(['view ',num2str(imgidx)])
    imgnum = imgnum+1;
    imgidx = imgidx+1;
end