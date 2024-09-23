clear all
close all

dataset = 1
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
params.ICL_DATA                 = 0;
params.multiK                   = 1;
params.syn                      = 0;

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

if(dataset == 1)
    ref_idx = 7;
    load('GenData\gamma1_ABC006_42n47.mat')
    Gamma1s_delta1 = load('GenData\gamma1_ABC006_42n47_delta1.mat').Gamma1s;
    params.HYPO1_VIEW_INDX          = 43;
    params.HYPO2_VIEW_INDX          = 48;
elseif(dataset == 2)
    ref_idx = 5;
    load('GenData\gamma1_ABC077_20n48.mat')
    params.HYPO1_VIEW_INDX          = 21;
    params.HYPO2_VIEW_INDX          = 49;
elseif(dataset == 3)
    ref_idx = 1;
    load('GenData\gamma1_ABC325_3n11.mat')
    params.HYPO1_VIEW_INDX          = 4;
    params.HYPO2_VIEW_INDX          = 12;
elseif(dataset == 4)
    ref_idx = 22;
    load('GenData\gamma1_ABC568_27n48.mat')
    params.HYPO1_VIEW_INDX          = 28;
    params.HYPO2_VIEW_INDX          = 49;
elseif(dataset == 5)
    ref_idx = 5;
    load('GenData\gamma1_ABC2211_41n47.mat')
    params.HYPO1_VIEW_INDX          = 42;
    params.HYPO2_VIEW_INDX          = 48;
end
recons_coor1 = Gamma1s';
recons_coor1_delta1 = Gamma1s_delta1';
abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
    T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_trans = R21afterBA'*(recons_coor1'-T21afterBA);
Gamma1s_trans_delta1 = R21afterBA'*(recons_coor1_delta1'-T21afterBA);

if(dataset == 1)
    load('GenData\gamma1_ABC006_22n26.mat')
    Gamma1s_delta1 = load('GenData\gamma1_ABC006_22n26_delta1.mat').Gamma1s;
    params.HYPO1_VIEW_INDX          = 23;
    params.HYPO2_VIEW_INDX          = 27;
elseif(dataset == 2)
    load('GenData\gamma1_ABC077_10n36.mat')
    params.HYPO1_VIEW_INDX          = 11;
    params.HYPO2_VIEW_INDX          = 37;
elseif(dataset == 3)
    load('GenData\gamma1_ABC325_5n48.mat')
    params.HYPO1_VIEW_INDX          = 6;
    params.HYPO2_VIEW_INDX          = 47;
elseif(dataset == 4)
    load('GenData\gamma1_ABC568_30n33.mat')
    params.HYPO1_VIEW_INDX          = 31;
    params.HYPO2_VIEW_INDX          = 33;
elseif(dataset == 5)
    load('GenData\gamma1_ABC2211_46n49.mat')
    params.HYPO1_VIEW_INDX          = 47;
    params.HYPO2_VIEW_INDX          = 50;
end
recons_coor2 = Gamma1s';
recons_coor2_delta1 = Gamma1s_delta1';
abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
    T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_trans1= R21afterBA'*(recons_coor2'-T21afterBA);
Gamma1s_trans1_delta1= R21afterBA'*(recons_coor2_delta1'-T21afterBA);

if(dataset == 1)
    load('GenData\gamma1_ABC006_6n8.mat')
    Gamma1s_delta1 = load('GenData\gamma1_ABC006_6n8_delta1.mat').Gamma1s;
    params.HYPO1_VIEW_INDX          = 7;
    params.HYPO2_VIEW_INDX          = 9;
elseif(dataset == 2)
    load('GenData\gamma1_ABC077_4n23.mat')
    params.HYPO1_VIEW_INDX          = 5;
    params.HYPO2_VIEW_INDX          = 24;
elseif(dataset == 3)
    load('GenData\gamma1_ABC325_0n16.mat')
    params.HYPO1_VIEW_INDX          = 1;
    params.HYPO2_VIEW_INDX          = 17;
elseif(dataset == 4)
    load('GenData\gamma1_ABC568_21n39.mat')
    params.HYPO1_VIEW_INDX          = 22;
    params.HYPO2_VIEW_INDX          = 40;
elseif(dataset == 5)
    load('GenData\gamma1_ABC2211_4n18.mat')
    params.HYPO1_VIEW_INDX          = 5;
    params.HYPO2_VIEW_INDX          = 19;
end
recons_coor3 = Gamma1s';
recons_coor3_delta1 = Gamma1s_delta1';

fignum = 1;
subnum = 1;
e3 = [0;0;1];
abs_R1afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C1afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);

recons_coorall = [recons_coor3; Gamma1s_trans'; Gamma1s_trans1'];
recons_coorall_delta1 = [recons_coor3_delta1; Gamma1s_trans_delta1'; Gamma1s_trans1_delta1'];
fig = 3;
sub = 1;
sbr = 2;
sbc = 5;
%{
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
    reprojection = PmatrixafterBA*[recons_coor1(:,1)';...
        recons_coor1(:,2)';...
        recons_coor1(:,3)';...
        ones(1,size(recons_coor1,1))];
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
    

% % % %     dist2Dall = dist2Dall';
% % % %     [dist2D1,dist2D_idx1]=mink(dist2Dall,1,2);
% % % %     edge_left  = find(dist2D1(:,1) >= 2*sqrt(2));
% % % % %     edge_left = dist2D_idx(fp_pt,:);
% % % %     leftedge{valiimg,1} = TO_Edges_1(edge_left,:);
% % % %     leftedgeidx{valiimg} = edge_left;
%     plot(TO_Edges_1(edge_left,1), TO_Edges_1(edge_left,2), 'y.', 'MarkerSize',4, 'LineWidth', 1);
    plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'y.', 'MarkerSize',4, 'LineWidth', 1);
    hold on;
    plot(x_coor(1,tp_pt), y_coor(1,tp_pt), 'g.', 'MarkerSize',3, 'LineWidth', 1);
    hold on;
    plot(x_coor(1,fp_pt), y_coor(1,fp_pt), 'r.', 'MarkerSize',3, 'LineWidth', 1);
%     legend ({'original edges' 'reprojected reconstruction edges tp' 'reprojected reconstruction edges fp'},...
%         'FontSize',12);
    title (['view ',num2str(valididx-1)])
    if(sub == 7 || sub == 10)
        legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
    end
    sub = sub +1;
    if(sub>(sbr*sbc))
        sub = 1;
            fig = fig;
    end
end
%}

Gamma_all = recons_coorall';
Gamma_all_delta1 = recons_coorall_delta1';
I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
    [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_trans1_1 = R21afterBA'*(recons_coorall'-T21afterBA);
Gamma1s_trans1_1_delta1 = R21afterBA'*(recons_coorall_delta1'-T21afterBA);

scale_min = [inf inf inf];
scale_max = [-inf -inf -inf];
if(dataset == 1)
    input_curves = load('MATData\curves0006.mat').curve_points;
    curve_MCS    = load('GenData\MCS_00000006.mat').recs;
elseif(dataset == 2)
    input_curves = load('MATData\curves0077.mat').curve_points;
    curve_MCS    = load('GenData\MCS_00000077.mat').recs;
elseif(dataset == 3)
    input_curves = load('MATData\curves0325.mat').curve_points;
    curve_MCS    = load('GenData\MCS_00000325.mat').recs;
    elseif(dataset == 4)
    input_curves = load('MATData\curves0568.mat').curve_points;
    curve_MCS    = load('GenData\MCS_00000568.mat').recs;
elseif(dataset == 5)
    input_curves = load('MATData\curves2211.mat').curve_points;
    curve_MCS    = load('GenData\MCS_00002211.mat').recs;
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
GT_all = [];
for i = 1:size(final_curves, 2)
    displacement = 0.5 - (point_location_max ./ 2);
    c = final_curves{i};
    c(:,1) = c(:,1) + displacement(1);
    c(:,2) = c(:,2) + displacement(2);
    c(:,3) = c(:,3) + displacement(3);
    final_curves{i} = c;
    GT_all          = [GT_all;c];
    p1 = plot3(final_curves{i}(:, 1), final_curves{i}(:, 2), final_curves{i}(:, 3), 'g', 'LineWidth', 1.5);
    hold on
end
hold off;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

hold on
MCS_all = [];
for i = 1:size(curve_MCS, 2)
    c_MCS = curve_MCS{i};
    MCS_all = [MCS_all; c_MCS];
    p2 = plot3(c_MCS(:, 1), c_MCS(:, 2), c_MCS(:, 3), 'b', 'LineWidth', 0.5);
    hold on
end

Gamma1s_trans1_1 = -Gamma1s_trans1_1;
Gamma1s_trans1_1_delta1 = -Gamma1s_trans1_1_delta1;
hold on
p3 = plot3(Gamma1s_trans1_1(1,:),       Gamma1s_trans1_1(2,:),       Gamma1s_trans1_1(3,:), ...
    'rx', 'MarkerSize',1, 'LineWidth', 1);
axis equal
title ('Visualization of ground truth 3D curve points, 3D curves from MCS and 3D reconstruction result','FontSize',12)
legend ([p1 p2 p3],{'3D curve points generated from ABC dataset(GT)', '3D curves from MCS', '3D recontruction result'},'FontSize',12)
%{
figure
subplot(2,2,1)
for i = 1:size(final_curves, 2)
    p1 = plot3(final_curves{i}(:, 1), final_curves{i}(:, 2), final_curves{i}(:, 3), 'g', 'LineWidth', 1.5);
    hold on
end

p3 = plot3(Gamma1s_trans1_1(1,:),       Gamma1s_trans1_1(2,:),       Gamma1s_trans1_1(3,:), ...
    'rx', 'MarkerSize',1, 'LineWidth', 1);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
title ('Visualization of ground truth 3D curve points and 3D reconstruction result','FontSize',12)

subplot(2,2,2)
for i = 1:size(final_curves, 2)
    p1 = plot3(final_curves{i}(:, 1), final_curves{i}(:, 2), final_curves{i}(:, 3), 'g', 'LineWidth', 1.5);
    hold on
end

for i = 1:size(curve_MCS, 2)
    c_MCS = curve_MCS{i};
    p2 = plot3(c_MCS(:, 1), c_MCS(:, 2), c_MCS(:, 3), 'b', 'LineWidth', 0.5);
    hold on
end
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
title ('Visualization of ground truth 3D curve points and 3D curves from MCS','FontSize',12)

subplot(2,2,3)
for i = 1:size(curve_MCS, 2)
    c_MCS = curve_MCS{i};
    p2 = plot3(c_MCS(:, 1), c_MCS(:, 2), c_MCS(:, 3), 'b', 'LineWidth', 0.5);
    hold on
end

p3 = plot3(Gamma1s_trans1_1(1,:),       Gamma1s_trans1_1(2,:),       Gamma1s_trans1_1(3,:), ...
    'rx', 'MarkerSize',1, 'LineWidth', 1);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
title ('Visualization of 3D curves from MCS and 3D reconstruction result','FontSize',12)
legend ([p1 p2 p3],{'3D curve points generated from ABC dataset(GT)', '3D curves from MCS', '3D recontruction result'},'FontSize',12)
%}

figure(1)
dataset =1
subplot(1,3,(dataset-1)*3+1)
for i = 1:size(curve_MCS, 2)
    c_MCS = curve_MCS{i};
    p2 = plot3(c_MCS(:, 1), c_MCS(:, 2), c_MCS(:, 3), 'b', 'LineWidth', 0.5);
    hold on
end
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
title ('3D Curves Sketch','FontSize',12)

subplot(1,3,(dataset-1)*3+2)
p3 = scatter3(Gamma1s_trans1_1(1,:), Gamma1s_trans1_1(2,:), Gamma1s_trans1_1(3,:), ...
        1,'MarkerEdgeColor','r','Marker','x');
% xlabel('x');
% ylabel('y');
% zlabel('z');
axis equal;
grid off
title ('3D Edge Sketch','FontSize',12)

subplot(1,3,3*dataset)
for i = 1:size(final_curves, 2)
    p1 = plot3(final_curves{i}(:, 1), final_curves{i}(:, 2), final_curves{i}(:, 3), 'g', 'LineWidth', 1.5);
    hold on
end
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
title ('Ground Truth','FontSize',12)


x_coorGT  = GT_all(:,1);
y_coorGT  = GT_all(:,2);
z_coorGT  = GT_all(:,3);

dx1 = ones(1,size(x_coorGT,1)).*Gamma1s_trans1_1(1,:)'-x_coorGT';
dy1 = ones(1,size(y_coorGT,1)).*Gamma1s_trans1_1(2,:)'-y_coorGT';
dz1 = ones(1,size(z_coorGT,1)).*Gamma1s_trans1_1(3,:)'-z_coorGT';

dist3Dall1 = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
[dist3D1,~]=mink(dist3Dall1,1,2);

dx = ones(1,size(x_coorGT,1)).*MCS_all(:,1)-x_coorGT';
dy = ones(1,size(y_coorGT,1)).*MCS_all(:,2)-y_coorGT';
dz = ones(1,size(z_coorGT,1)).*MCS_all(:,3)-z_coorGT';

dist3Dall = sqrt(dx.^2 + dy.^2 + dz.^2);
[dist3D,~]=mink(dist3Dall,1,2);

figure
subplot(1,2,1)
histogram(dist3D1,'BinWidth', 0.005,"Normalization","probability")
xlabel ('distance to the GT 3D vurve point', 'FontSize',10)
ylabel ('Number of  Occurrence', 'FontSize',10)
xlim([0,max([dist3D1;dist3D])])
ylim([0,1])
title ('Distribution of minimum distance between 3D edge sketch and GT', 'FontSize',10)

subplot(1,2,2)
histogram(dist3D,'BinWidth', 0.005,"Normalization","probability")
xlabel ('distance to the GT 3D vurve point', 'FontSize',10)
ylabel ('Number of  Occurrence', 'FontSize',10)
title ('Distribution of minimum distance between 3D curve sketch and GT', 'FontSize',10)
xlim([0,max([dist3D1;dist3D])])
ylim([0,1])


figure
subplot(1,2,1)
p3 = plot3(Gamma1s_trans1_1(1,:),       Gamma1s_trans1_1(2,:),       Gamma1s_trans1_1(3,:), ...
    'bx', 'MarkerSize',1, 'LineWidth', 1);
axis equal
xlim([min(Gamma1s_trans1_1_delta1(1,:)) max(Gamma1s_trans1_1_delta1(1,:))])
ylim([min(Gamma1s_trans1_1_delta1(2,:)) max(Gamma1s_trans1_1_delta1(2,:))])
zlim([min(Gamma1s_trans1_1_delta1(3,:)) max(Gamma1s_trans1_1_delta1(3,:))])
title ('3D Edge Sketch','FontSize',12)
subplot(1,2,2)
p3 = plot3(Gamma1s_trans1_1_delta1(1,:),       Gamma1s_trans1_1_delta1(2,:),       Gamma1s_trans1_1_delta1(3,:), ...
    'rx', 'MarkerSize',1, 'LineWidth', 1);
axis equal
xlim([min(Gamma1s_trans1_1_delta1(1,:)) max(Gamma1s_trans1_1_delta1(1,:))])
ylim([min(Gamma1s_trans1_1_delta1(2,:)) max(Gamma1s_trans1_1_delta1(2,:))])
zlim([min(Gamma1s_trans1_1_delta1(3,:)) max(Gamma1s_trans1_1_delta1(3,:))])
title ('3D Edge Sketch with delta=0.1','FontSize',12)

Edge_all = Gamma1s_trans1_1;
Edge_all_delta1 = Gamma1s_trans1_1_delta1;

% 
% for(dataset = 2:5)
% subplot(5,3,(dataset-1)*3+1)
% title ''
% subplot(5,3,(dataset-1)*3+2)
% title ''
% subplot(5,3,3*dataset)
% title ''
% end
% legend ([p1 p2 p3],{'3D curve points generated from ABC dataset(GT)', '3D curves from MCS', '3D recontruction result'},'FontSize',12)