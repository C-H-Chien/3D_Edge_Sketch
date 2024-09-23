% Epipolar wedges method is used in this main code.
clear all; 
fileID = fopen('GenData/pairededge_ABC0006_6n8_t32to0excludehypo1n2_delta03_theta15_N4.txt','r');
% fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
B = reshape(A,50,size(A,1)/50);
B=B+1;
paired_edge = B';
% 
% load('MATData/TO_Edges_10_sigma15.mat');
% load('MATData/imageArray_10.mat')
% load('MATData/R_matrix.mat')
% load('MATData/T_matrix.mat')
% % load('MATData/R_matrix_10.mat')
% % load('MATData/T_matrix_10.mat')
% load('MATData/K_10.mat')
% thinner_pair = load('GenData/wedge0077.mat').paired_edge;

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

% load('MATData/TO_Edges_ABC00000077.mat');
% load('MATData/imageArray_ABC00000077.mat')
% load('MATData/R_matrix_ABC00000077.mat')
% load('MATData/T_matrix_ABC00000077.mat')
% load('MATData/K_ABC00000077.mat')
% collect_third_order_edges = collect_third_order_edges';
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
        TO_Edgels_Orientations = collect_third_order_edges{i,6}(:,3);

        Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
        for j = 1:size(TO_Edgels_Orientations, 1)
            Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
                sin(TO_Edgels_Orientations(j,1))];
        end
        TO_Edges_Tangents{i,1} = Tgts;
    end

else
    for i = 1:params.NUM_OF_IMGS
        TO_Edges_Tangents{i,1} = collect_third_order_edges{i,6}(:,3:4);
    end
end

% tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
% tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};

%> Array of view indices
view_idx = [];
% figure;
for i = 1:params.NUM_OF_IMGS; view_idx = [view_idx; i]; end
validation_view_indices = view_idx;
validation_view_indices([params.HYPO1_VIEW_INDX; ...
                         params.HYPO2_VIEW_INDX],:) = [];

%> Fetch edgels of the two hypothesis views
% load('F:\T-less\edges50_1.mat');
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,6};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,6};
Gamma1s         = [];
GammaTangent_3D = [];

hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
rows      = size(hypo_img2, 1);
cols      = size(hypo_img2, 2);

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

finalEdgePair    = [];
% paired_edge      = [];

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

for(i = 1 : size(paired_edge,1))
    edgel_HYPO1 = TO_Edges_HYPO1(paired_edge(i,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge(i,2), :);
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    [edgels_HYPO2_corrected, ...
     edgels_HYPO1_corrected, e_i] = edgelsHYPO2correct(Epipolar_Coeffs, ...
                                                       edgel_HYPO2, ...
                                                       edgel_HYPO1, ...
                                                       R_matrix, ...
                                                       T_matrix, ...
                                                       params, K);
    edgels_HYPO2_final = edgels_HYPO2_corrected;
    edgels_HYPO1_final = edgels_HYPO1_corrected;
    if(params.multiK == 1)
                    pt1 = invK1 * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK2 * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                else
                    pt1 = invK * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
    end
    pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
    %pts_meters = [edgel_HYPO1(1, 1), edgels_HYPO2(finalPairIndx, 1); edgel_HYPO1(1, 2), edgels_HYPO2(finalPairIndx, 2)];
    Gamma = linearTriangulation(2, pts_meters, R21, T21);
    GammaTgt_3D     = get3DOrientation(edgels_HYPO1_final,edgels_HYPO2_final, K, params, R21);
    Gamma1s = [Gamma1s, Gamma];
    GammaTangent_3D = [GammaTangent_3D, GammaTgt_3D];
    Gamma = Gamma*1000;
%     hold on;
%     plot3(Gamma(1,:), Gamma(2,:), Gamma(3,:), 'b.');
%     hold on;
%     plot3([Gamma(1,1)-GammaTgt_3D(1,1) Gamma(1,1)+GammaTgt_3D(1,1)], ...
%           [Gamma(2,1)-GammaTgt_3D(2,1) Gamma(2,1)+GammaTgt_3D(2,1)], ...
%           [Gamma(3,1)-GammaTgt_3D(3,1) Gamma(3,1)+GammaTgt_3D(3,1)], 'b');
end

I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_orin = -R21afterBA'*(Gamma1s-T21afterBA);
Gamma1s_ore_orin1 = -R21afterBA'*(GammaTangent_3D);

figure;
plot3(Gamma1s_orin(1,:)*100, Gamma1s_orin(2,:)*100, Gamma1s_orin(3,:)*100, 'b.');
axis equal;
title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels, threshold starts from ', num2str(32),')']) 


% figure;
% for(id = 1 : size(paired_edge,1))
% hold on;
% colorcur = rand(1,3);
% quiver3(Gamma1s_orin(1,id)*500, Gamma1s_orin(2,id)*500, Gamma1s_orin(3,id)*500,Gamma1s_ore_orin(1,id), Gamma1s_ore_orin(2,id), Gamma1s_ore_orin(3,id),...
%     "AutoScale","off",'ShowArrowHead','off','Color',colorcur)
% end
% axis equal
% view(3)

figure;
for(id = 1 : size(paired_edge,1))
hold on;
colorcur = rand(1,3);
plot3(Gamma1s_orin(1,id)*500, Gamma1s_orin(2,id)*500, Gamma1s_orin(3,id)*500, ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',7)
Gamma1s_ore_orin = Gamma1s_ore_orin1/3;
plot3([Gamma1s_orin(1,id)*500-Gamma1s_ore_orin(1,id), Gamma1s_orin(1,id)*500+Gamma1s_ore_orin(1,id)], ...
      [Gamma1s_orin(2,id)*500-Gamma1s_ore_orin(2,id), Gamma1s_orin(2,id)*500+Gamma1s_ore_orin(2,id)], ...
      [Gamma1s_orin(3,id)*500-Gamma1s_ore_orin(3,id), Gamma1s_orin(3,id)*500+Gamma1s_ore_orin(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
axis equal
view(3)