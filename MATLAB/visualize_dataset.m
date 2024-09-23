clear all
close all
load('MATData/TO_Edges_ABC00000006.mat');
load('MATData/imageArray_ABC00000006.mat')

% load('F:\231017\TO_Edges_rgbd_dataset_freiburg3_cabinet_validation.mat');
% load('F:\231017\imageAparams.rray_rgbd_dataset_freiburg3_cabinet_validation.mat')

params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 12;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*5;
params.ICL_DATA                 = 0;
params.multiK                   = 1;
params.cols                     = size(double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX})), 2);
vi = 1;

% i = 3;
% len = 1;
% figure
% TO_Edges_HYPO1 = collect_third_order_edges{i,size(collect_third_order_edges,2)};
% imshow(imageArray{i});
% hold on;
% plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',7);
% Tgts = zeros(size(TO_Edges_HYPO1, 1), 2);
% for j = 1:size(TO_Edges_HYPO1, 1)
%     edgel_tgt1 = [cos(TO_Edges_HYPO1(j,3)), ...
%         sin(TO_Edges_HYPO1(j,3))];
%     x1 = TO_Edges_HYPO1(j,1)+len*edgel_tgt1(1,1);
%     x2 = TO_Edges_HYPO1(j,1)-len*edgel_tgt1(1,1);
%     y1 = TO_Edges_HYPO1(j,2)+len*edgel_tgt1(1,2);
%     y2 = TO_Edges_HYPO1(j,2)-len*edgel_tgt1(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'c', 'LineWidth', 0.25);
% end

i = 7;

TO_Edges_HYPO1 = collect_third_order_edges{i,size(collect_third_order_edges,2)};
imshow(imageArray{i});
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'm.', 'MarkerSize',1.5, 'LineWidth', 2);



for(i = 1:25)
    figure(16)
    TO_Edges_HYPO1 = collect_third_order_edges{i,size(collect_third_order_edges,2)};
    subplot(5,5,i)
    imshow(imageArray{i});
    hold on;
    plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'm.', 'MarkerSize',1.5, 'LineWidth', 2);
    title(['view ',num2str(i-1)])
    figure(17)
    TO_Edges_HYPO1 = collect_third_order_edges{i+25,size(collect_third_order_edges,2)};
    subplot(5,5,i)
    imshow(imageArray{i+25});
    hold on;
    plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'm.', 'MarkerSize',1.5, 'LineWidth', 2);
    title(['view ', num2str(i+25-1)])
end

% figure
% for(i = 1:10)
%     TO_Edges_HYPO1 = collect_third_order_edges{i,size(collect_third_order_edges,2)};
%     subplot(3,4,i)
%     imshow(imageArray{i});
%     hold on;
%     plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',0.25);
%     title(['view ', num2str(i)])
% end