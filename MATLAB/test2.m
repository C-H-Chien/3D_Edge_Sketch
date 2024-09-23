% close all
clear all

load('MATData/TO_Edges_10.mat');
load('MATData/imageArray_10.mat')
load('MATData/R_matrix_10.mat')
load('MATData/T_matrix_10.mat')
load('MATData/K_10.mat')
load('MATData/imageName_10.mat')

load('GenData/paired_tless_6n16cppThinWedge.mat')

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

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

%> Dataset and sequence paths
datasetFolder = 'F:\T-less\';
sequenceName  = '10\';
%> Data paths and file names
FilePath_Depths   = strcat(datasetFolder, sequenceName, 'depth\');

%> Fetch edgels of the two hypothesis views
TO_Edges_HYPO1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
TO_Edges_HYPO2 = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};

%> Read Images
HYPO1_IMG = imageArray{params.HYPO1_VIEW_INDX};
HYPO2_IMG = imageArray{params.HYPO2_VIEW_INDX};

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
tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};

%> Get Relative Pose
[R21, T21, E21, F21] = getRelativePose(R_matrix, T_matrix, params, K);

% %> Display
% figure;
% subplot(1,2,1)
% imshow(rgb2gray(HYPO1_IMG));
% hold on;
% plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',4, 'LineWidth', 1);
% title 'hypothesis view 1'
% 
% subplot(1,2,2)
% imshow(rgb2gray(HYPO2_IMG)); 
% hold on;
% plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',4, 'LineWidth', 1);
% title 'hypothesis view 2'
% hold on;

correct_num  = 0;
correct_pair = [];
wrong_num    = 0;
wrong_pair   = [];
all_dist     = [];
wrong_pair_to_investigate = [];
wrong_pair40 = [];
wrong_pair20 = [];
wrong_pair10 = [];
% % % load('MATData/ambi_edge4133_99.mat')
for (idx = 1000 : size(paired_edge,1))
% % %     if(isempty(find(ambi_edge == paired_edge(idx,1))) ~= 1)
% % %         continue;
% % %     end
%     if(isempty(find(ambi_edge == paired_edge(idx,1))) ~= 1)
%         continue
%     end
%     paired_edge(idx,1)
    edgel_HYPO1 = round(TO_Edges_HYPO1(paired_edge(idx,1), 1:2));
    edgel_HYPO1_1 = TO_Edges_HYPO1(paired_edge(idx,1), 1:2);
    edgel_tgt1  = tgt1(paired_edge(idx,1), 1:2);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge(idx,2), 1:2);
    edgel_tgt2  = tgt2(paired_edge(idx,2), 1:2);

    %> Compute the coefficients of the epipolar line
    coeffs = F21 * [edgel_HYPO1_1'; 1];
    Epipolar_Coeffs.A = coeffs(1,1);
    Epipolar_Coeffs.B = coeffs(2,1);
    Epipolar_Coeffs.C = coeffs(3,1);

    %> Read Depths
    HYPO1_Depth_dir = strcat(FilePath_Depths, imageName{params.HYPO1_VIEW_INDX}, '.png');
    HYPO2_Depth_dir = strcat(FilePath_Depths, imageName{params.HYPO2_VIEW_INDX}, '.png');
    HYPO1_Depth = double(imread(HYPO1_Depth_dir));
    HYPO2_Depth = double(imread(HYPO2_Depth_dir));
%     HYPO1_Depth = HYPO1_Depth ./ 5000;
%     HYPO2_Depth = HYPO2_Depth ./ 5000;
    HYPO1_Depth = HYPO1_Depth .*0.1;
    HYPO2_Depth = HYPO2_Depth .*0.1;

    %> Get ground truth correspondence point gamma2 on the second image from
    %  gamma1 on the first image
    e3 = [0;0;1];
    if(params.multiK == 1)
        gamma1 = invK1 * [edgel_HYPO1'; 1];
    else
        gamma1 = invK * [edgel_HYPO1'; 1];
    end
    rho1   = HYPO1_Depth(edgel_HYPO1(2), edgel_HYPO1(1));
    gamma2 = (rho1*R21*gamma1 + T21)/((e3'*R21*gamma1)*rho1 + e3'*T21);
    if(params.multiK == 1)
        pt2    = K2*gamma2;
    else
        pt2    = K*gamma2;
    end

    pt2_2d = pt2(1:2,1)';

    distance = norm(edgel_HYPO2-pt2_2d);
    all_dist = [all_dist; distance];
    len = 3;

%     realpt2idxraw = find(TO_Edges_HYPO2(:,1) >= floor(pt2(1)) & TO_Edges_HYPO2(:,1) <= ceil(pt2(1)) & ...
%          TO_Edges_HYPO2(:,2) >= floor(pt2(2)) & TO_Edges_HYPO2(:,2) <= ceil(pt2(2)));
%     if(isempty(realpt2idxraw))
%         continue;
%     elseif(size(realpt2idxraw,1)>1)
% % %         dist_pt2 = [];
% % %         for(ipt2 = 1: size(TO_Edges_HYPO2,1))
% % %             dist_pt2 = [dist_pt2; norm(TO_Edges_HYPO2(ipt2, 1:2)-pt2_2d)];
% % %         end
% % %         [minidist,pt2idx] = mink(dist_pt2,1);
% % %         realpt2idx = pt2idx;
%         if(minidist >= 2)
%             continue;
%         end
%     else
%         realpt2idx = realpt2idxraw;
%     end

    %
    %> Display
    figure(1)
    subplot(1,2,1)
    imshow(rgb2gray(HYPO1_IMG));
    hold on;
    plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',4, 'LineWidth', 1);
    hold on;
    plot(edgel_HYPO1_1(1), edgel_HYPO1_1(2), 'mo', 'MarkerSize',10, 'LineWidth', 1);
    plot(edgel_HYPO1_1(1), edgel_HYPO1_1(2), 'mx', 'MarkerSize',10, 'LineWidth', 1);
    x1 = edgel_HYPO1_1(1,1)+len*edgel_tgt1(1,1);
    x2 = edgel_HYPO1_1(1,1)-len*edgel_tgt1(1,1);
    y1 = edgel_HYPO1_1(1,2)+len*edgel_tgt1(1,2);
    y2 = edgel_HYPO1_1(1,2)-len*edgel_tgt1(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'm', 'LineWidth', 1);

    subplot(1,2,2)
    imshow(rgb2gray(HYPO2_IMG));
    hold on;
    plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',4, 'LineWidth', 1);
    hold on;
%     pt2 = TO_Edges_HYPO2(realpt2idx, 1:2);
    cols = size(HYPO2_IMG, 2);
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
    line([1, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 0.5); hold on;
    line([pt2(1), edgel_HYPO2(1)], [pt2(2), edgel_HYPO2(2)], 'Color', 'w', 'LineWidth', 1);
    plot(pt2(1), pt2(2), 'go', 'MarkerSize',10, 'LineWidth', 1);
    plot(pt2(1), pt2(2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
%     edgel_tgt21 = tgt2(realpt2idx, 1:2);
%     x1 = pt2(1)+len*edgel_tgt21(1,1);
%     x2 = pt2(1)-len*edgel_tgt21(1,1);
%     y1 = pt2(2)+len*edgel_tgt21(1,2);
%     y2 = pt2(2)-len*edgel_tgt21(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'g', 'LineWidth', 1);
    plot(edgel_HYPO2(1), edgel_HYPO2(2), 'ro', 'MarkerSize',10, 'LineWidth', 1);
    plot(edgel_HYPO2(1), edgel_HYPO2(2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
    x1 = edgel_HYPO2(1,1)+len*edgel_tgt2(1,1);
    x2 = edgel_HYPO2(1,1)-len*edgel_tgt2(1,1);
    y1 = edgel_HYPO2(1,2)+len*edgel_tgt2(1,2);
    y2 = edgel_HYPO2(1,2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'r', 'LineWidth', 1);

    hold on;
    set(gcf,'color','w');
    %}
% % 
% %     if(distance <=2)
% %         correct_num  = correct_num+1;
% %         correct_pair = [correct_pair; paired_edge(idx,:)];
% %     else
% %         wrong_num  = wrong_num+1;
% %         wrong_pair = [wrong_pair; paired_edge(idx,:)];
% %         if(distance >= 40)
% %             wrong_pair40 = [wrong_pair40; paired_edge(idx,:)];
% %         end
% %         if(distance >= 20)
% %             wrong_pair20 = [wrong_pair20; paired_edge(idx,:)];
% %         end
% %         if(distance >= 10)
% %             wrong_pair10 = [wrong_pair10; paired_edge(idx,:)];
% %         end
% %     end
end
figure;
histogram(all_dist, 'BinWidth',1,'Normalization','probability')
ylabel 'probability'
xlabel 'error'
title 'distribution of euclidean distance'