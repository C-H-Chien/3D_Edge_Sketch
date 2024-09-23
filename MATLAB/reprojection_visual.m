dataset = 1;
if(dataset == 1)
    load('MATData/TO_Edges_ABC00000006.mat');
    load('MATData/imageArray_ABC00000006.mat')
    load('MATData/R_matrix_ABC00000006.mat')
    load('MATData/T_matrix_ABC00000006.mat')
    load('MATData/K_ABC00000006.mat')
elseif(dataset == 2)
    load('MATData/TO_Edges_ABC00000077.mat');
    load('MATData/imageArray_ABC00000077.mat')
    load('MATData/R_matrix_ABC00000077.mat')
    load('MATData/T_matrix_ABC00000077.mat')
    load('MATData/K_ABC00000077.mat')
elseif(dataset == 3)
    load('MATData/TO_Edges_ABC00000325.mat');
    load('MATData/imageArray_ABC00000325.mat')
    load('MATData/R_matrix_ABC00000325.mat')
    load('MATData/T_matrix_ABC00000325.mat')
    load('MATData/K_ABC00000325.mat')
elseif(dataset == 4)
    load('MATData/TO_Edges_ABC00000568.mat');
    load('MATData/imageArray_ABC00000568.mat')
    load('MATData/R_matrix_ABC00000568.mat')
    load('MATData/T_matrix_ABC00000568.mat')
    load('MATData/K_ABC00000568.mat')
elseif(dataset == 5)
    load('MATData/TO_Edges_ABC00002211.mat');
    load('MATData/imageArray_ABC00002211.mat')
    load('MATData/R_matrix_ABC00002211.mat')
    load('MATData/T_matrix_ABC00002211.mat')
    load('MATData/K_ABC00002211.mat')
end
% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 49;
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


leftedge = [];
leftedgeidx = [];

TO_Edges_Tangents = cell(params.NUM_OF_IMGS, 1);
for i = 1:params.NUM_OF_IMGS
    TO_Edgels_Orientations = collect_third_order_edges{i,size(collect_third_order_edges,2)}(:,3);

    Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
    for j = 1:size(TO_Edgels_Orientations, 1)
        Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
            sin(TO_Edgels_Orientations(j,1))];
    end
    TO_Edges_Tangents{i,1} = Tgts;
end
view_idx = [];
for i = 1:params.NUM_OF_IMGS; view_idx = [view_idx; i]; end
validation_view_indices = view_idx;
validation_view_indices([params.HYPO1_VIEW_INDX; ...
    params.HYPO2_VIEW_INDX],:) = [];
if(dataset == 1)
ref_idx = 7;
load('GenData\gamma1_ABC006_42n47.mat')
params.HYPO1_VIEW_INDX          = 43;
params.HYPO2_VIEW_INDX          = 48;
recons_coor1 = Gamma1s';
elseif(dataset == 2)
ref_idx = 5;
load('GenData\gamma1_ABC077_20n48.mat')
params.HYPO1_VIEW_INDX          = 21;
params.HYPO2_VIEW_INDX          = 49;
recons_coor1 = Gamma1s';
elseif(dataset == 3)
ref_idx = 1;
load('GenData\gamma1_ABC325_3n11.mat')
params.HYPO1_VIEW_INDX          = 4;
params.HYPO2_VIEW_INDX          = 12;
recons_coor1 = Gamma1s';
elseif(dataset == 4)
ref_idx = 22;
load('GenData\gamma1_ABC568_27n48.mat')
params.HYPO1_VIEW_INDX          = 28;
params.HYPO2_VIEW_INDX          = 49;
recons_coor1 = Gamma1s';
elseif(dataset == 5)
ref_idx = 5;
load('GenData\gamma1_ABC2211_41n47.mat')
params.HYPO1_VIEW_INDX          = 42;
params.HYPO2_VIEW_INDX          = 48;
recons_coor1 = Gamma1s';
end
% load("GenData\pair_ore15_delta03_6n12_sigma15.mat")

abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
                 T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_trans = R21afterBA'*(recons_coor1'-T21afterBA);

if(dataset == 1)
load('GenData\gamma1_ABC006_22n26.mat')
recons_coor2 = Gamma1s';
params.HYPO1_VIEW_INDX          = 23;
params.HYPO2_VIEW_INDX          = 27;
elseif(dataset == 2)
load('GenData\gamma1_ABC077_10n36.mat')
recons_coor2 = Gamma1s';
params.HYPO1_VIEW_INDX          = 11;
params.HYPO2_VIEW_INDX          = 37;
elseif(dataset == 3)
load('GenData\gamma1_ABC325_5n48.mat')
recons_coor2 = Gamma1s';
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 47;
elseif(dataset == 4)
load('GenData\gamma1_ABC568_30n33.mat')
recons_coor2 = Gamma1s';
params.HYPO1_VIEW_INDX          = 31;
params.HYPO2_VIEW_INDX          = 33;
elseif(dataset == 5)
load('GenData\gamma1_ABC2211_46n49.mat')
recons_coor2 = Gamma1s';
params.HYPO1_VIEW_INDX          = 47;
params.HYPO2_VIEW_INDX          = 50;
end
% load("GenData\pair_ore15_delta03_6n12_sigma15.mat")

abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
                 T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_trans1= R21afterBA'*(recons_coor2'-T21afterBA);

% % % % % % % load('GenData\gamma1_ABC32t5_3n11.mat')
% % % % % % % figure(12)
% % % % % % % subplot(1,3,1)
% % % % % % % plot3(Gamma1s(1,:),       -Gamma1s(2,:),       -Gamma1s(3,:), ...
% % % % % % %       'b.', 'MarkerSize',3, 'LineWidth', 1);
% % % % % % % axis equal;
% % % % % % % title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels)',...
% % % % % % %         ' hypo view: view 6&8 multiple round']) 
% % % % % % % 
% % % % % % % subplot(1,3,2)
% % % % % % % plot3(Gamma1s_trans(1,:), -Gamma1s_trans(2,:), -Gamma1s_trans(3,:), ...
% % % % % % %       'm.', 'MarkerSize',3, 'LineWidth', 1);
% % % % % % % axis equal;
% % % % % % % title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels)',...
% % % % % % %         ' hypo view: view 42&47 single round']) 
% % % % % % % 
% % % % % % % subplot(1,3,3)
% % % % % % % plot3(Gamma1s_trans1(1,:), -Gamma1s_trans1(2,:), -Gamma1s_trans1(3,:), ...
% % % % % % %       'm.', 'MarkerSize',3, 'LineWidth', 1);
% % % % % % % axis equal;
% % % % % % % title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels)',...
% % % % % % %         ' hypo view: view 22&26 single round']) 

% idx6_12 = [];
% idx9_15 = [];
% dist_all = [];
% for(pt3d_i = 1:size(Gamma1s,2))
%     curr_3dpt6n12 = Gamma1s(:,pt3d_i);
%     dxdydz        = Gamma1s_trans - curr_3dpt6n12;
%     dxdydzsqr     = dxdydz.^2;
%     distsqr       = sum(dxdydzsqr);
%     [min_d,min_i] = mink(distsqr,1);
%     dist_all      = [dist_all; min_d];
%     if(min_d<=0.003)
%         idx6_12 = [idx6_12;pt3d_i];
%         idx9_15 = [idx9_15;min_i];
%     end
% end
% 
% v6_12only = 1:size(Gamma1s,2);
% v6_12only(idx6_12) = [];
% v9_15only = 1:size(Gamma1s_trans,2);
% v9_15only(idx9_15) = [];
% 
% 
% figure(11)
% subplot(1,2,1)
% plot3([0], [0], [0], ...
%       'b.', 'MarkerSize',10, 'LineWidth', 3);
% hold on;
% plot3([0] , [0], [0], ...
%       'm.', 'MarkerSize',10, 'LineWidth', 3);
% hold on;
% plot3([0] , [0], [0], ...
%       'g.', 'MarkerSize',10, 'LineWidth', 3);
% 
% plot3(Gamma1s(1,v6_12only),       -Gamma1s(2,v6_12only), ...
%       -Gamma1s(3,v6_12only),       'b.', 'MarkerSize',3, 'LineWidth', 1);
% hold on;
% plot3(Gamma1s(1,idx6_12),         -Gamma1s(2,idx6_12), ...
%       -Gamma1s(3,idx6_12),         'g.', 'MarkerSize',5, 'LineWidth', 2);
% hold on;
% plot3(Gamma1s_trans(1,v9_15only), -Gamma1s_trans(2,v9_15only), ...
%       -Gamma1s_trans(3,v9_15only), 'm.', 'MarkerSize',3, 'LineWidth', 1);
% hold on;
% plot3(Gamma1s_trans(1,idx9_15),   -Gamma1s_trans(2,idx9_15), ...
%       -Gamma1s_trans(3,idx9_15),   'g.', 'MarkerSize',5, 'LineWidth', 2);
% axis equal;
% xlim([min([Gamma1s(1,:),Gamma1s_trans(1,:)]),   ...
%       max([Gamma1s(1,:),Gamma1s_trans(1,:)])])
% ylim([min([-Gamma1s(2,:),-Gamma1s_trans(2,:)]), ...
%       max([-Gamma1s(2,:),-Gamma1s_trans(2,:)])])
% zlim([min([-Gamma1s(3,:),-Gamma1s_trans(3,:)]), ...
%       max([-Gamma1s(3,:),-Gamma1s_trans(3,:)])])
% legend 'only in multi' 'only in  single' 'in both'
% title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels, ',...
%         'combination of two results)']) 
% 
% 
% subplot(1,2,2)
% plot3(Gamma1s(1,idx6_12),         -Gamma1s(2,idx6_12), ...
%       -Gamma1s(3,idx6_12),         'g.', 'MarkerSize',5, 'LineWidth', 2);
% hold on;
% plot3(Gamma1s_trans(1,idx9_15),   -Gamma1s_trans(2,idx9_15), ...
%       -Gamma1s_trans(3,idx9_15),   'g.', 'MarkerSize',5, 'LineWidth', 2);
% axis equal;
% xlim([min([Gamma1s(1,:),Gamma1s_trans(1,:)]),   ...
%       max([Gamma1s(1,:),Gamma1s_trans(1,:)])])
% ylim([min([-Gamma1s(2,:),-Gamma1s_trans(2,:)]), ...
%       max([-Gamma1s(2,:),-Gamma1s_trans(2,:)])])
% zlim([min([-Gamma1s(3,:),-Gamma1s_trans(3,:)]), ...
%       max([-Gamma1s(3,:),-Gamma1s_trans(3,:)])])
% title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels,',...
%         ' in both two results)']) 
% 

if(dataset == 1)
load('GenData\gamma1_ABC006_6n8.mat')
recons_coor1 = Gamma1s';
params.HYPO1_VIEW_INDX          = 7;
params.HYPO2_VIEW_INDX          = 9;
elseif(dataset == 2)
load('GenData\gamma1_ABC077_4n23.mat')
recons_coor1 = Gamma1s';
params.HYPO1_VIEW_INDX          = 5;
params.HYPO2_VIEW_INDX          = 24;
elseif(dataset == 3)
load('GenData\gamma1_ABC325_0n16.mat')
recons_coor1 = Gamma1s';
params.HYPO1_VIEW_INDX          = 1;
params.HYPO2_VIEW_INDX          = 17;
elseif(dataset == 4)
load('GenData\gamma1_ABC568_21n39.mat')
recons_coor1 = Gamma1s';
params.HYPO1_VIEW_INDX          = 22;
params.HYPO2_VIEW_INDX          = 40;
elseif(dataset == 5)
load('GenData\gamma1_ABC2211_4n18.mat')
recons_coor1 = Gamma1s';
params.HYPO1_VIEW_INDX          = 5;
params.HYPO2_VIEW_INDX          = 19;
end
fignum = 1;
subnum = 1;
e3 = [0;0;1];
abs_R1afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C1afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
%{
figure(1);
subplot(1,2,1)
TO_Edges_1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,size(collect_third_order_edges,2)};
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'm.', 'MarkerSize',5, 'LineWidth', 5);
hold on
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'g.', 'MarkerSize',5, 'LineWidth', 5);
hold on
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'rx', 'MarkerSize',5, 'LineWidth', 1);
imshow(imageArray{params.HYPO1_VIEW_INDX});
hold on;
plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'm.', 'MarkerSize',3, 'LineWidth', 1);
Kmatrix       = K(:,:,params.HYPO1_VIEW_INDX);
extrafterBA           = zeros(3,4);
extrafterBA(1:3, 1:3) = [1,0,0;0,1,0;0,0,1];
extrafterBA(1:3, 4)   = [0;0;0];
PmatrixafterBA        = Kmatrix*extrafterBA;
reprojection = PmatrixafterBA*[recons_coor1(:,1)';...
                               recons_coor1(:,2)';...
                               recons_coor1(:,3)';...
                               ones(1,size(recons_coor1,1))];
x_coor       = reprojection(1,:)./reprojection(3,:);
y_coor       = reprojection(2,:)./reprojection(3,:);
dx = ones(1,size(TO_Edges_1,1)).*x_coor' - TO_Edges_1(:,1)';
dy = ones(1,size(TO_Edges_1,1)).*y_coor' - TO_Edges_1(:,2)';
dist2Dall = sqrt(dx.^2 + dy.^2);
dist2D = min(dist2Dall, [], 2);
fp_pt1  = find(dist2D(:,1) >= 2);

tp_pt1 = 1:size(x_coor,2);
tp_pt1(fp_pt1) = [];
hold on;
plot(x_coor(1,tp_pt1), y_coor(1,tp_pt1), 'g.', 'MarkerSize',4, 'LineWidth', 1);
hold on;
plot(x_coor(1,fp_pt1), y_coor(1,fp_pt1), 'rx', 'MarkerSize',6, 'LineWidth', 0.5);
% legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
title 'hypothesis 1 view 7'

subplot(1,2,2);
TO_Edges_1 = collect_third_order_edges{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'm.', 'MarkerSize',5, 'LineWidth', 5);
hold on
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'g.', 'MarkerSize',5, 'LineWidth', 5);
hold on
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'rx', 'MarkerSize',5, 'LineWidth', 1);
imshow(imageArray{params.HYPO2_VIEW_INDX});
hold on;
% TO_Edges_1 = collect_third_order_edges{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'm.', 'MarkerSize',3, 'LineWidth', 1);
abs_R2afterBA =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
    T_matrix(:,params.HYPO2_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Kmatrix       = K(:,:,params.HYPO2_VIEW_INDX);
extrafterBA           = zeros(3,4);
extrafterBA(1:3, 1:3) = R21afterBA;
extrafterBA(1:3, 4)   = T21afterBA;
PmatrixafterBA        = Kmatrix*extrafterBA;
reprojection = PmatrixafterBA*[recons_coor1(:,1)';...
                               recons_coor1(:,2)';...
                               recons_coor1(:,3)';...
                               ones(1,size(recons_coor1,1))];
x_coor       = reprojection(1,:)./reprojection(3,:);
y_coor       = reprojection(2,:)./reprojection(3,:);
dx = ones(1,size(TO_Edges_1,1)).*x_coor' - TO_Edges_1(:,1)';
dy = ones(1,size(TO_Edges_1,1)).*y_coor' - TO_Edges_1(:,2)';
dist2Dall = sqrt(dx.^2 + dy.^2);
dist2D = min(dist2Dall, [], 2);
fp_pt  = find(dist2D(:,1) >= 2);

tp_pt = 1:size(x_coor,2);
tp_pt(fp_pt) = [];
hold on;
plot(x_coor(1,tp_pt), y_coor(1,tp_pt), 'g.', 'MarkerSize',4, 'LineWidth', 1);
hold on;
plot(x_coor(1,fp_pt), y_coor(1,fp_pt), 'rx', 'MarkerSize',6, 'LineWidth', 0.5);
% legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
title 'hypothesis 2 view 9'
%}

%{
err_pt = fp_pt(13,1);

figure(2);
subplot(1,2,1)
TO_Edges_1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,size(collect_third_order_edges,2)};
paired_edge_fp = paired_edge(err_pt,1);
imshow(imageArray{params.HYPO1_VIEW_INDX});
hold on;
plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'c.', 'MarkerSize',3, 'LineWidth', 1);
hold on;
plot(TO_Edges_1(paired_edge_fp{1},1), TO_Edges_1(paired_edge_fp{1},2), 'rx', 'MarkerSize',8, 'LineWidth', 2);
% legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
title 'hypothesis 1 view 6'

subplot(1,2,2);
TO_Edges_1 = collect_third_order_edges{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
paired_edge_fp = paired_edge(err_pt,2);
imshow(imageArray{params.HYPO2_VIEW_INDX});
hold on;
plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'c.', 'MarkerSize',3, 'LineWidth', 1);
hold on;
plot(TO_Edges_1(paired_edge_fp{1},1), TO_Edges_1(paired_edge_fp{1},2), 'rx', 'MarkerSize',8, 'LineWidth', 2);
title 'hypothesis 2 view 12'

figure(1);
subplot(1,2,1)
TO_Edges_1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,size(collect_third_order_edges,2)};
plot(TO_Edges_1(1,1), TO_Edges_1(1,2), 'm.', 'MarkerSize',5, 'LineWidth', 5);
Kmatrix       = K(:,:,params.HYPO1_VIEW_INDX);
extrafterBA           = zeros(3,4);
extrafterBA(1:3, 1:3) = [1,0,0;0,1,0;0,0,1];
extrafterBA(1:3, 4)   = [0;0;0];
PmatrixafterBA        = Kmatrix*extrafterBA;
reprojection = PmatrixafterBA*[recons_coor1(:,1)';...
                               recons_coor1(:,2)';...
                               recons_coor1(:,3)';...
                               ones(1,size(recons_coor1,1))];
x_coor       = reprojection(1,:)./reprojection(3,:);
y_coor       = reprojection(2,:)./reprojection(3,:);
dx = ones(1,size(TO_Edges_1,1)).*x_coor' - TO_Edges_1(:,1)';
dy = ones(1,size(TO_Edges_1,1)).*y_coor' - TO_Edges_1(:,2)';
dist2Dall = sqrt(dx.^2 + dy.^2);
dist2D = min(dist2Dall, [], 2);
fp_pt1  = find(dist2D(:,1) >= 2);

tp_pt1 = 1:size(x_coor,2);
tp_pt1(fp_pt1) = [];
hold on;
plot(x_coor(1,tp_pt1), y_coor(1,tp_pt1), 'g.', 'MarkerSize',4, 'LineWidth', 1);
hold on;
plot(x_coor(1,fp_pt1), y_coor(1,fp_pt1), 'rx', 'MarkerSize',6, 'LineWidth', 0.5);
hold on;
% plot(x_coor(1,err_pt), y_coor(1,err_pt), 'ro', 'MarkerSize',8, 'LineWidth', 1);
% legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
title 'hypothesis 1 view 6'

subplot(1,2,2);
TO_Edges_1 = collect_third_order_edges{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
imshow(imageArray{params.HYPO2_VIEW_INDX});
hold on;
% TO_Edges_1 = collect_third_order_edges{params.HYPO2_VIEW_INDX,size(collect_third_order_edges,2)};
plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'm.', 'MarkerSize',3, 'LineWidth', 1);
abs_R2afterBA =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
    T_matrix(:,params.HYPO2_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Kmatrix       = K(:,:,params.HYPO2_VIEW_INDX);
extrafterBA           = zeros(3,4);
extrafterBA(1:3, 1:3) = R21afterBA;
extrafterBA(1:3, 4)   = T21afterBA;
PmatrixafterBA        = Kmatrix*extrafterBA;
reprojection = PmatrixafterBA*[recons_coor1(:,1)';...
                               recons_coor1(:,2)';...
                               recons_coor1(:,3)';...
                               ones(1,size(recons_coor1,1))];
x_coor       = reprojection(1,:)./reprojection(3,:);
y_coor       = reprojection(2,:)./reprojection(3,:);
dx = ones(1,size(TO_Edges_1,1)).*x_coor' - TO_Edges_1(:,1)';
dy = ones(1,size(TO_Edges_1,1)).*y_coor' - TO_Edges_1(:,2)';
dist2Dall = sqrt(dx.^2 + dy.^2);
dist2D = min(dist2Dall, [], 2);
fp_pt  = find(dist2D(:,1) >= 2);

tp_pt = 1:size(x_coor,2);
tp_pt(fp_pt) = [];
hold on;
plot(x_coor(1,tp_pt), y_coor(1,tp_pt), 'g.', 'MarkerSize',4, 'LineWidth', 1);
hold on;
plot(x_coor(1,fp_pt), y_coor(1,fp_pt), 'rx', 'MarkerSize',6, 'LineWidth', 0.5);
hold on;
plot(x_coor(1,err_pt), y_coor(1,err_pt), 'ro', 'MarkerSize',8, 'LineWidth', 1);
% sgtitle (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels, orientation different ≤ ', num2str(params.SUPPORT_OREN_THRESH),'°)'])
%}
recons_coor1_all = [recons_coor1; Gamma1s_trans'; Gamma1s_trans1'];
fig = 3;
sub = 1;
sbr = 2;
sbc = 5;
%
for(valiimg = 1:50)
    figure(fig)
    subplot(sbr,sbc,sub)
%     valididx = validation_view_indices(valiimg,1);
valididx = valiimg;
    imshow(imageArray{valididx});
    hold on;
    TO_Edges_1 = collect_third_order_edges{valididx,size(collect_third_order_edges,2)};
    
    abs_R2afterBA =  R_matrix(:,:,valididx);
    abs_C2afterBA = -R_matrix(:,:,valididx)' *...
        T_matrix(:,valididx);
    R21afterBA    = abs_R2afterBA * abs_R1afterBA';
    T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
    Kmatrix       = K(:,:,valididx);
    extrafterBA           = zeros(3,4);
    extrafterBA(1:3, 1:3) = R21afterBA;
    extrafterBA(1:3, 4)   = T21afterBA;
    PmatrixafterBA        = Kmatrix*extrafterBA;
    reprojection = PmatrixafterBA*[recons_coor1_all(:,1)';...
        recons_coor1_all(:,2)';...
        recons_coor1_all(:,3)';...
        ones(1,size(recons_coor1_all,1))];
    x_coor       = reprojection(1,:)./reprojection(3,:);
    y_coor       = reprojection(2,:)./reprojection(3,:);

    dx = ones(1,size(TO_Edges_1,1)).*x_coor' - TO_Edges_1(:,1)';
    dy = ones(1,size(TO_Edges_1,1)).*y_coor' - TO_Edges_1(:,2)';
    dist2Dall = sqrt(dx.^2 + dy.^2);
    [dist2D,dist2D_idx]=mink(dist2Dall,1,2);
    fp_pt  = find(dist2D(:,1) >= 2*sqrt(2));
%     edge_left = dist2D_idx(fp_pt,:);
%     leftedge{valiimg} = TO_Edges_1(edge_left,:);
%     leftedgeidx{valiimg} = edge_left;

    tp_pt = 1:size(x_coor,2);
    tp_pt(fp_pt) = [];
    edge_ignore = dist2D_idx(tp_pt,:);
    dist2Dall = sqrt(dx.^2 + dy.^2);

    [dist2D_1,dist2D_idx_1]=mink(dist2Dall',1,2);
    fp_pt_1  = find(dist2D_1(:,1) >= 2*sqrt(2));
%     edge_left = dist2D_idx(fp_pt,:);
%     leftedge{valiimg} = TO_Edges_1(edge_left,:);
%     leftedgeidx{valiimg} = edge_left;

    tp_pt_1 = 1:size(dist2Dall,2);
    tp_pt_1(fp_pt_1) = [];
%     edge_ignore = dist2D_idx(tp_pt,:);
    

% % % %     dist2Dall = dist2Dall';
% % % %     [dist2D1,dist2D_idx1]=mink(dist2Dall,1,2);
% % % %     edge_left  = find(dist2D1(:,1) >= 2*sqrt(2));
% % % % %     edge_left = dist2D_idx(fp_pt,:);
% % % %     leftedge{valiimg,1} = TO_Edges_1(edge_left,:);
% % % %     leftedgeidx{valiimg} = edge_left;
%     plot(TO_Edges_1(edge_left,1), TO_Edges_1(edge_left,2), 'y.', 'MarkerSize',4, 'LineWidth', 1);
    plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'y.', 'MarkerSize',4, 'LineWidth', 1);
    plot(TO_Edges_1(fp_pt_1,1), TO_Edges_1(fp_pt_1,2), 'm.', 'MarkerSize',4, 'LineWidth', 1);
%     hold on;
%     plot(x_coor(1,tp_pt), y_coor(1,tp_pt), 'g.', 'MarkerSize',3, 'LineWidth', 1);
%     hold on;
%     plot(x_coor(1,fp_pt), y_coor(1,fp_pt), 'r.', 'MarkerSize',3, 'LineWidth', 1);
%     legend ({'original edges' 'reprojected reconstruction edges tp' 'reprojected reconstruction edges fp'},...
%         'FontSize',12);
    title (['view ',num2str(valididx-1)])
%     if(sub == 7 || sub == 10)
%         legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
%     end
    sub = sub +1;
    if(sub>(sbr*sbc))
        sub = 1;
            fig = fig+1;
    end
end
%}

% load('GenData\gamma1_ore15_delta03_9n15_sigma15_left.mat')

Gamma_all = recons_coor1_all';
%
I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma_all_trans1 = -R21afterBA'*(recons_coor1'-T21afterBA);
Gamma_all_trans2 = -R21afterBA'*(Gamma1s_trans-T21afterBA);
Gamma_all_trans3 = -R21afterBA'*(Gamma1s_trans1-T21afterBA);

% figure;
% plot3(Gamma_all_trans1(1,:),       -Gamma_all_trans1(2,:),       -Gamma_all_trans1(3,:), ...
%       'b.', 'MarkerSize',3, 'LineWidth', 1);
% axis equal;
% title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels)',...
%         ' accumulative']) 

figure
plot3(Gamma_all_trans1(1,:), Gamma_all_trans1(2,:), Gamma_all_trans1(3,:), ...
      'b.', 'MarkerSize',3, 'LineWidth', 1);
hold on;
plot3(Gamma_all_trans2(1,:), Gamma_all_trans2(2,:), Gamma_all_trans2(3,:), ...
      'm.', 'MarkerSize',3, 'LineWidth', 1);
hold on;
plot3(Gamma_all_trans3(1,:), Gamma_all_trans3(2,:), Gamma_all_trans3(3,:), ...
      'g.', 'MarkerSize',3, 'LineWidth', 1);
axis equal;
title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels)',...
        ' accumulative']) 
%}

I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_trans1_1 = R21afterBA'*(recons_coor1_all'-T21afterBA);
%{
figure;
plot3(Gamma1s_trans1_1(1,:),       -Gamma1s_trans1_1(2,:),       -Gamma1s_trans1_1(3,:), ...
      'b.', 'MarkerSize',3, 'LineWidth', 1);
axis equal;
title (['3D reconstruction result (Δ = ', num2str(params.delta),'pixels)',...
        ' accumulative']) 
%}

scale_min = [inf inf inf];
scale_max = [-inf -inf -inf];
if(dataset == 1)
input_curves = load('MATData\curves0006.mat').curve_points;
elseif(dataset == 2)
input_curves = load('MATData\curves0077.mat').curve_points;
elseif(dataset == 3)
input_curves = load('MATData\curves0325.mat').curve_points;
elseif(dataset == 4)
input_curves = load('MATData\curves0568.mat').curve_points;
elseif(dataset == 5)
input_curves = load('MATData\curves2211.mat').curve_points;
end

for i = 1:size(input_curves, 2)
    c = input_curves{i};
    c_min = min(c);
    c_max = max(c);
    scale_min = min([scale_min; c_min]);
    scale_max = max([scale_max; c_max]);
end
%> calculate the scale factor. Using this factor to contain the object in a
%1 by 1 by 1 bounding box
factor = max([scale_max - scale_min]');
factor = 1/factor;

%> Move all curve points to the first quatriple and scale up/down so that 
%  one of the axis is in the interval [0,1]
point_location_max = [-inf -inf -inf];
shift_curves = input_curves;
for i = 1:size(shift_curves, 2)
    c = shift_curves{i};
    c(:,1) = (c(:,1) - scale_min(1))*factor;
    c(:,2) = (c(:,2) - scale_min(2))*factor;
    c(:,3) = (c(:,3) - scale_min(3))*factor;
    shift_curves{i} = c;

    c_max = max(c);
    point_location_max = max([point_location_max; c_max]);
end

%> Shift the entire curve points to center at (0.5, 0.5)
figure;
final_curves = shift_curves;
for i = 1:size(final_curves, 2)
    displacement = 0.5 - (point_location_max ./ 2);
    c = final_curves{i};
    c(:,1) = c(:,1) + displacement(1);
    c(:,2) = c(:,2) + displacement(2);
    c(:,3) = c(:,3) + displacement(3);
    final_curves{i} = c;

    plot3(c(:, 1), c(:, 2), c(:, 3));
    hold on
end
hold off;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

Gamma1s_trans1_1 = -Gamma1s_trans1_1;

% scale_min_1 = [inf inf inf];
% scale_max_1 = [-inf -inf -inf];
% scale_min_1 = min(Gamma1s_trans1_1');
% scale_max_1 = max(Gamma1s_trans1_1');
% movept      =  abs(scale_max_1-scale_min_1)/2;
% Gamma1s_trans1_2(1,:) = Gamma1s_trans1_1(1,:) - scale_min_1(1,1);
% Gamma1s_trans1_2(2,:) = Gamma1s_trans1_1(2,:) - movept(1,2);
% Gamma1s_trans1_2(3,:) = Gamma1s_trans1_1(3,:) - scale_min_1(1,3);
hold on
plot3(Gamma1s_trans1_1(1,:),       Gamma1s_trans1_1(2,:),       Gamma1s_trans1_1(3,:), ...
'b.', 'MarkerSize',1, 'LineWidth', 1);
axis equal