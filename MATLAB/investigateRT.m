% Epipolar wedges method is used in this main code.
clear all; 
% close all;

% load all files needed
% load('MATData/TO_Edges_10_sigma15.mat');
% load('MATData/imageArray_10.mat')
% load('MATData/R_matrix.mat')
% load('MATData/T_matrix.mat')
% load('MATData/K_10.mat')

load('MATData/TO_Edges_ABC00000006.mat');
load('MATData/imageArray_ABC00000006.mat')
load('MATData/R_matrix_ABC00000006.mat')
load('MATData/T_matrix_ABC00000006.mat')
load('MATData/K_ABC00000006.mat')

% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 7;
params.HYPO2_VIEW_INDX          = 9;
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
params.circle                   = 55;
params.parallelangle            = 15;
% params.SUPPORT_OREN_THRESH_DEG  = 15;

hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));

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
if(params.syn == 0)
    for i = 1:params.NUM_OF_IMGS
        TO_Edgels_Orientations = collect_third_order_edges{i,size(collect_third_order_edges,2)}(:,3);

        Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
        for j = 1:size(TO_Edgels_Orientations, 1)
            Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
                sin(TO_Edgels_Orientations(j,1))];
        end
        TO_Edges_Tangents{i,size(collect_third_order_edges,2)} = Tgts;
    end

else
    for i = 1:params.NUM_OF_IMGS
        TO_Edges_Tangents{i,size(collect_third_order_edges,2)} = collect_third_order_edges{i,size(collect_third_order_edges,2)}(:,3:4);
    end
end

tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,size(collect_third_order_edges,2)};
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};

%> Array of view indices
view_idx = [];
for i = 1:params.NUM_OF_IMGS; view_idx = [view_idx; i]; end
validation_view_indices = view_idx;
validation_view_indices([params.HYPO1_VIEW_INDX; ...
                         params.HYPO2_VIEW_INDX],:) = [];

%> Fetch edgels of the two hypothesis views
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,size(collect_third_order_edges,2)};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
Gamma1s         = [];
GammaTangent_3D = [];

hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
rows      = size(hypo_img2, 1);
cols      = size(hypo_img2, 2);

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

finalEdgePair    = [];
paired_edge      = [];

% get the orientation list between hypo1 and hypo2
[ore_list1bar_all, ...
 ore_list2_sorted, ...
 epipole_pix_view1, ...
 epipole_pix_view2]= getOreList_New(TO_Edges_HYPO1, TO_Edges_HYPO2, ...
                                     R_matrix, T_matrix, ...
                                     params.HYPO1_VIEW_INDX, ...
                                     params.HYPO2_VIEW_INDX, params, K);

% figure(1)
% subplot(1,2,1);
% imshow(uint8(hypo_img1));
% hold on;
% plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'y.', 'MarkerSize',5, 'LineWidth', 2);
% title 'hypothesis view 1'
% subplot(1,2,2);
% imshow(uint8(hypo_img2));
% hold on;
% plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'y.', 'MarkerSize',5, 'LineWidth', 2);
% title 'hypothesis view 2'

% pipeline starts
cnt_sup   = 0;
sup_num   = [];
ambi_edge = [];
each_round= [];

edge_idx = 555;
edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
Apixel = coeffs(1,1);
Bpixel = coeffs(2,1);
Cpixel = coeffs(3,1);
Epipolar_Coeffs.A = Apixel;
Epipolar_Coeffs.B = Bpixel;
Epipolar_Coeffs.C = Cpixel;
%
figure(2)
tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,size(collect_third_order_edges,2)};
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
subplot(1,2,1);
imshow(uint8(hypo_img1));
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
dxdy_hypo1 = edgel_HYPO1(1, 1:2) - epipole_pix_view1(1:2,1)';
slope1     = dxdy_hypo1(1,2)/dxdy_hypo1(1,1);
yintersect = epipole_pix_view1(2,1)-slope1*epipole_pix_view1(1,1);
yMin       = yintersect;
yMax       = slope1*cols+yintersect;
hold on;
line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 1);
title 'hypothesis view 1'

subplot(1,2,2);
imshow(uint8(hypo_img2));
title 'hypothesis view 2'
hold on;
plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 1);
T_x = @(T)[0,      -T(3,1),  T(2,1); ...
           T(3,1),  0,      -T(1,1); ...
          -T(2,1),  T(1,1),  0];

eulZYX = rotm2eul(R21);
percentageR = 0.975:0.005:1.025;
percentageT = 0.975:0.005:1.025;
for(preturbi = 1:size(percentageT,2))
    eulZYX_pret = [eulZYX(1,1)*1, ...
                   eulZYX(1,2)*percentageR(1,preturbi), ...
                   eulZYX(1,3)*1];
%     eulZYX_pret = eulZYX;
%     T21_pre     = [T21(1,1)*1; ...
%                    T21(2,1)*percentageT(1,preturbi); ...
%                    T21(3,1)*1];
    T21_pre     = T21;
    if(percentageR(1,preturbi) == 1)
        continue;
    end
    R21_pret    = eul2rotm(eulZYX_pret);
    E   = T_x(T21_pre) * R21_pret;
    % Calculate fundamental matrix
    if(params.multiK == 1)
        K1 = K(:,:,params.HYPO1_VIEW_INDX);
        K2 = K(:,:,params.HYPO2_VIEW_INDX);
        invK1 = inv(K1);
        invK2 = inv(K2);
        Fpre   = invK2'*E*invK1;
    else
        invK = inv(K);
        Fpre   = invK'*E*invK;
    end
    coeffs1 = Fpre * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs1(1,1);
    Bpixel = coeffs1(2,1);
    Cpixel = coeffs1(3,1);
    Epipolar_Coeffs1.A = Apixel;
    Epipolar_Coeffs1.B = Bpixel;
    Epipolar_Coeffs1.C = Cpixel;
    yMin = -Epipolar_Coeffs1.C./Epipolar_Coeffs1.B;
    yMax = (-Epipolar_Coeffs1.C - Epipolar_Coeffs1.A*cols) ./ Epipolar_Coeffs1.B;
    line([0, cols], [yMin, yMax], 'Color', 'r', 'LineWidth', 0.5);
end
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 1);
